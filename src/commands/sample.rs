use anyhow::{bail, Context, Result};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;

use crate::{
    cli::SampleArgs,
    formats::{SeqFormat, SeqRecord},
    io, utils,
};

pub fn run(args: SampleArgs) -> Result<()> {
    if args.num.is_none() && args.rate.is_none() {
        bail!("one of --num or --rate is required");
    }
    if let Some(rate) = args.rate {
        if !(0.0..=1.0).contains(&rate) {
            bail!("--rate must be in [0,1]");
        }
    }
    let mut rng = ChaCha20Rng::seed_from_u64(args.seed);

    let paired_in1 = args.input1.as_deref().or(args.in1.as_deref());
    let paired_in2 = args.input2.as_deref().or(args.in2.as_deref());
    utils::validate_input_mode("sample", args.io.input.as_deref(), paired_in1, paired_in2)?;
    if let (Some(in1), Some(in2)) = (paired_in1, paired_in2) {
        log_info(&args, "sample paired-end mode");
        let out2 = args
            .output2
            .as_deref()
            .context("paired sample requires --output2")?;
        if let Some(rate) = args.rate {
            log_info(&args, "sample by rate");
            let outputted = sample_paired_rate_streaming(&args, in1, in2, out2, rate, &mut rng)?;
            log_info(&args, &format!("{outputted} read pairs outputted"));
            return Ok(());
        }
        log_info(&args, "sample by number");
        let n = args.num.unwrap();
        let (o1, o2, processed) = sample_paired_num_streaming(
            in1,
            in2,
            &args.io.compression,
            n,
            &mut rng,
            args.allow_unpaired,
        )?;
        io::write_record_pair_parallel(
            &args.io.output,
            out2,
            SeqFormat::Fastq,
            &args.io.compression,
            &o1,
            &o2,
        )?;
        log_info(
            &args,
            &format!(
                "{processed} read pairs processed, {} read pairs outputted",
                o1.len()
            ),
        );
        return Ok(());
    }

    run_single_streaming(args, &mut rng)
}

fn run_single_streaming(args: SampleArgs, rng: &mut ChaCha20Rng) -> Result<()> {
    let pb = make_progress(args.progress);
    let in_path = args.io.input.as_deref();
    let (fmt, reader) = io::open_seq_reader(in_path, &args.io.format, &args.io.compression)?;

    if let Some(rate) = args.rate {
        log_info(&args, "sample by rate");
        return sample_rate_streaming(&args, reader, fmt, rate, rng, &pb);
    }

    log_info(&args, "sample by number");
    let n = args.num.unwrap();
    let (out, processed) = reservoir_from_reader(reader, fmt, n, rng, &pb)?;
    finish_progress(&pb);
    log_info(
        &args,
        &format!(
            "{processed} sequences processed, {} sequences outputted",
            out.len()
        ),
    );
    io::write_records(&args.io.output, fmt, &args.io.compression, &out)
}

fn sample_paired_rate_streaming(
    args: &SampleArgs,
    in1: &str,
    in2: &str,
    out2: &str,
    rate: f64,
    rng: &mut ChaCha20Rng,
) -> Result<usize> {
    let w1 = io::open_writer(&args.io.output)?;
    let w1 = io::wrap_compress_for_streams(w1, &args.io.output, &args.io.compression, 2)?;
    let w2 = io::open_writer(out2)?;
    let w2 = io::wrap_compress_for_streams(w2, out2, &args.io.compression, 2)?;
    let mut w1 = io::buffered_writer(w1);
    let mut w2 = io::buffered_writer(w2);

    let mut outputted = 0usize;
    let mut invalid_preview = Vec::new();
    let mut invalid_count = 0usize;
    io::for_each_fastq_pair(in1, in2, &args.io.compression, |pair| {
        match pair {
            io::FastqPair::Both(a, b) => {
                let a_id = a.id_str()?;
                let b_id = b.id_str()?;
                if pair_key(a_id) != pair_key(b_id) {
                    invalid_count += 1;
                    if invalid_preview.len() < 10 {
                        invalid_preview.push(format!("R1={a_id} R2={b_id}"));
                    }
                    return Ok(true);
                }
                if rng.random::<f64>() <= rate {
                    io::write_record_ref(&mut w1, SeqFormat::Fastq, &a)?;
                    io::write_record_ref(&mut w2, SeqFormat::Fastq, &b)?;
                    outputted += 1;
                }
            }
            io::FastqPair::Left(a) => {
                invalid_count += 1;
                if invalid_preview.len() < 10 {
                    invalid_preview.push(format!("unmatched in R1: {}", a.id_str()?));
                }
            }
            io::FastqPair::Right(b) => {
                invalid_count += 1;
                if invalid_preview.len() < 10 {
                    invalid_preview.push(format!("unmatched in R2: {}", b.id_str()?));
                }
            }
        }
        Ok(true)
    })?;
    report_invalid_pairs(invalid_count, invalid_preview, args.allow_unpaired)?;
    Ok(outputted)
}

fn sample_paired_num_streaming(
    in1: &str,
    in2: &str,
    compression: &crate::cli::CompressionArg,
    n: usize,
    rng: &mut ChaCha20Rng,
    allow_unpaired: bool,
) -> Result<(Vec<SeqRecord>, Vec<SeqRecord>, usize)> {
    if n == 0 {
        return Ok((Vec::new(), Vec::new(), 0));
    }
    let mut output = Vec::with_capacity(n);
    let mut valid_seen = 0usize;
    let mut invalid_preview = Vec::new();
    let mut invalid_count = 0usize;
    let mut w = 0.0f64;
    let mut gap_remaining = 0usize;
    io::for_each_fastq_pair(in1, in2, compression, |pair| {
        match pair {
            io::FastqPair::Both(a, b) => {
                let a_id = a.id_str()?;
                let b_id = b.id_str()?;
                if pair_key(a_id) != pair_key(b_id) {
                    invalid_count += 1;
                    if invalid_preview.len() < 10 {
                        invalid_preview.push(format!("R1={a_id} R2={b_id}"));
                    }
                    return Ok(true);
                }
                valid_seen += 1;
                if output.len() < n {
                    output.push((a.to_owned_record()?, b.to_owned_record()?));
                    if output.len() == n {
                        w = (sample_open01(rng).ln() / n as f64).exp();
                        gap_remaining = sample_gap(w, rng);
                    }
                } else if gap_remaining > 0 {
                    gap_remaining -= 1;
                } else {
                    let j = rng.random_range(0..n);
                    output[j] = (a.to_owned_record()?, b.to_owned_record()?);
                    w *= (sample_open01(rng).ln() / n as f64).exp();
                    gap_remaining = sample_gap(w, rng);
                }
            }
            io::FastqPair::Left(a) => {
                invalid_count += 1;
                if invalid_preview.len() < 10 {
                    invalid_preview.push(format!("unmatched in R1: {}", a.id_str()?));
                }
            }
            io::FastqPair::Right(b) => {
                invalid_count += 1;
                if invalid_preview.len() < 10 {
                    invalid_preview.push(format!("unmatched in R2: {}", b.id_str()?));
                }
            }
        }
        Ok(true)
    })?;
    report_invalid_pairs(invalid_count, invalid_preview, allow_unpaired)?;
    let mut out1 = Vec::with_capacity(output.len());
    let mut out2 = Vec::with_capacity(output.len());
    for (a, b) in output {
        out1.push(a);
        out2.push(b);
    }
    Ok((out1, out2, valid_seen))
}

fn report_invalid_pairs(count: usize, preview: Vec<String>, allow_unpaired: bool) -> Result<()> {
    if count == 0 {
        return Ok(());
    }
    let preview = preview.join("; ");
    if !allow_unpaired {
        bail!(
            "paired inputs contain {count} invalid/unpaired records; examples: {preview}. rerun with --allow-unpaired to continue"
        );
    }
    eprintln!(
        "warning: paired inputs contain {count} invalid/unpaired records; examples: {preview}. unpaired records will be skipped"
    );
    Ok(())
}

fn pair_key(id: &str) -> &str {
    id.strip_suffix("/1")
        .or_else(|| id.strip_suffix("/2"))
        .unwrap_or(id)
}

fn sample_rate_streaming(
    args: &SampleArgs,
    reader: io::DynBufReader,
    fmt: SeqFormat,
    rate: f64,
    rng: &mut ChaCha20Rng,
    pb: &Option<ProgressBar>,
) -> Result<()> {
    let writer = io::open_writer(&args.io.output)?;
    let writer = io::wrap_compress(writer, &args.io.output, &args.io.compression)?;
    let writer = io::buffered_writer(writer);
    let mut processed = 0usize;
    let mut outputted = 0usize;
    let mut writer = writer;
    io::for_each_record_from_reader(reader, fmt, |rec| {
        processed += 1;
        if rng.random::<f64>() <= rate {
            io::write_record_ref(&mut writer, fmt, &rec)?;
            outputted += 1;
        }
        update_progress(pb);
        Ok(())
    })?;
    finish_progress(pb);
    if !args.io.quiet {
        eprintln!("[INFO] {processed} sequences processed, {outputted} sequences outputted");
    }
    Ok(())
}

fn reservoir_from_reader(
    reader: io::DynBufReader,
    fmt: SeqFormat,
    n: usize,
    rng: &mut ChaCha20Rng,
    pb: &Option<ProgressBar>,
) -> Result<(Vec<SeqRecord>, usize)> {
    if n == 0 {
        return Ok((Vec::new(), 0));
    }
    let mut out = Vec::with_capacity(n);
    let mut processed = 0usize;
    let mut w = 0.0f64;
    let mut gap_remaining = 0usize;

    io::for_each_record_from_reader(reader, fmt, |rec| {
        processed += 1;
        if out.len() < n {
            out.push(rec.to_owned_record()?);
            update_progress(pb);
            if out.len() == n {
                w = (sample_open01(rng).ln() / n as f64).exp();
                gap_remaining = sample_gap(w, rng);
            }
            return Ok(());
        }

        if gap_remaining > 0 {
            gap_remaining -= 1;
            update_progress(pb);
            return Ok(());
        }

        update_progress(pb);

        let j = rng.random_range(0..n);
        out[j] = rec.to_owned_record()?;

        w *= (sample_open01(rng).ln() / n as f64).exp();
        gap_remaining = sample_gap(w, rng);
        Ok(())
    })?;
    Ok((out, processed))
}

fn sample_gap(w: f64, rng: &mut ChaCha20Rng) -> usize {
    (sample_open01(rng).ln() / (1.0 - w).ln()).floor() as usize
}

fn sample_open01(rng: &mut ChaCha20Rng) -> f64 {
    loop {
        let u = rng.random::<f64>();
        if (0.0..1.0).contains(&u) {
            return u;
        }
    }
}

fn log_info(args: &SampleArgs, msg: &str) {
    if !args.io.quiet {
        eprintln!("[INFO] {msg}");
    }
}

fn make_progress(enabled: bool) -> Option<ProgressBar> {
    if !enabled {
        return None;
    }
    let pb = ProgressBar::with_draw_target(None, ProgressDrawTarget::stderr());
    pb.set_style(
        ProgressStyle::with_template("{spinner:.green} {pos} records processed")
            .expect("valid progress template"),
    );
    pb.enable_steady_tick(std::time::Duration::from_millis(120));
    Some(pb)
}

fn update_progress(pb: &Option<ProgressBar>) {
    if let Some(pb) = pb {
        pb.inc(1);
    }
}

fn finish_progress(pb: &Option<ProgressBar>) {
    if let Some(pb) = pb {
        pb.finish_and_clear();
    }
}
