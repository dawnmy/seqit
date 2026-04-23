use anyhow::{bail, Context, Result};
use bio::io::{fasta, fastq};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;
use std::io::BufReader;

use crate::{
    cli::{FormatArg, SampleArgs},
    formats::{SeqFormat, SeqRecord},
    io, pairs,
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
    if let (Some(in1), Some(in2)) = (paired_in1, paired_in2) {
        log_info(&args, "sample paired-end mode");
        let r1 = io::read_records(Some(in1), SeqFormat::Fastq, &args.io.compression)?;
        let r2 = io::read_records(Some(in2), SeqFormat::Fastq, &args.io.compression)?;
        let (r1, r2) = pairs::prepare_paired_records(r1, r2, args.allow_unpaired)?;
        let out2 = args
            .output2
            .as_deref()
            .context("paired sample requires --output2")?;
        if let Some(rate) = args.rate {
            log_info(&args, "sample by rate");
            let mut o1 = Vec::new();
            let mut o2 = Vec::new();
            for (a, b) in r1.into_iter().zip(r2.into_iter()) {
                if rng.random::<f64>() <= rate {
                    o1.push(a);
                    o2.push(b);
                }
            }
            io::write_records(&args.io.output, SeqFormat::Fastq, &args.io.compression, &o1)?;
            io::write_records(out2, SeqFormat::Fastq, &args.io.compression, &o2)?;
            log_info(&args, &format!("{} read pairs outputted", o1.len()));
            return Ok(());
        }
        log_info(&args, "sample by number");
        let n = args.num.unwrap();
        let mut idx: Vec<usize> = (0..r1.len()).collect();
        use rand::seq::SliceRandom;
        idx.shuffle(&mut rng);
        idx.truncate(n.min(idx.len()));
        idx.sort_unstable();
        let mut o1 = Vec::with_capacity(idx.len());
        let mut o2 = Vec::with_capacity(idx.len());
        for i in idx {
            o1.push(r1[i].clone());
            o2.push(r2[i].clone());
        }
        io::write_records(&args.io.output, SeqFormat::Fastq, &args.io.compression, &o1)?;
        io::write_records(out2, SeqFormat::Fastq, &args.io.compression, &o2)?;
        log_info(&args, &format!("{} read pairs outputted", o1.len()));
        return Ok(());
    }

    run_single_streaming(args, &mut rng)
}

fn run_single_streaming(args: SampleArgs, rng: &mut ChaCha20Rng) -> Result<()> {
    let pb = make_progress(args.progress);
    let in_path = args.io.input.as_deref();
    let fmt = if matches!(args.io.format, FormatArg::Auto) && matches!(in_path, None | Some("-")) {
        // stdin auto-detection must consume once.
        log_info(
            &args,
            "loading all sequences into memory for stdin auto-detect...",
        );
        let (fmt, recs) =
            io::read_records_with_format(in_path, &args.io.format, &args.io.compression)?;
        let out = if let Some(rate) = args.rate {
            log_info(&args, "sample by rate");
            recs.into_iter()
                .filter(|_| rng.random::<f64>() <= rate)
                .collect()
        } else {
            log_info(&args, "sample by number");
            reservoir_from_iter(recs.into_iter(), args.num.unwrap(), rng)
        };
        finish_progress(&pb);
        log_info(&args, &format!("{} sequences outputted", out.len()));
        return io::write_records(&args.io.output, fmt, &args.io.compression, &out);
    } else {
        io::resolve_seq_format(in_path, &args.io.format, &args.io.compression)?
    };

    if let Some(rate) = args.rate {
        log_info(&args, "sample by rate");
        return sample_rate_streaming(
            in_path,
            fmt,
            rate,
            &args.io.output,
            &args.io.compression,
            rng,
            &pb,
            args.io.quiet,
        );
    }

    log_info(&args, "sample by number");
    let n = args.num.unwrap();
    let (out, processed) = match fmt {
        SeqFormat::Fasta => {
            let r = io::open_reader(in_path)?;
            let r = io::wrap_decompress(r, in_path, &args.io.compression)?;
            let fa = fasta::Reader::new(BufReader::new(r));
            let iter = fa.records().map(|r| {
                let r = r?;
                Ok(SeqRecord {
                    id: r.id().to_string(),
                    desc: r.desc().map(|d| d.to_string()),
                    seq: r.seq().to_vec(),
                    qual: None,
                })
            });
            reservoir_from_result_iter(iter, n, rng, &pb)?
        }
        SeqFormat::Fastq => {
            let r = io::open_reader(in_path)?;
            let r = io::wrap_decompress(r, in_path, &args.io.compression)?;
            let fq = fastq::Reader::new(BufReader::new(r));
            let iter = fq.records().map(|r| {
                let r = r?;
                Ok(SeqRecord {
                    id: r.id().to_string(),
                    desc: r.desc().map(|d| d.to_string()),
                    seq: r.seq().to_vec(),
                    qual: Some(r.qual().to_vec()),
                })
            });
            reservoir_from_result_iter(iter, n, rng, &pb)?
        }
        _ => bail!("sample currently supports FASTA/FASTQ input"),
    };
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

fn sample_rate_streaming(
    in_path: Option<&str>,
    fmt: SeqFormat,
    rate: f64,
    output: &str,
    compression: &crate::cli::CompressionArg,
    rng: &mut ChaCha20Rng,
    pb: &Option<ProgressBar>,
    quiet: bool,
) -> Result<()> {
    let writer = io::open_writer(output)?;
    let writer = io::wrap_compress(writer, output, compression)?;
    let mut processed = 0usize;
    let mut outputted = 0usize;
    match fmt {
        SeqFormat::Fasta => {
            let r = io::open_reader(in_path)?;
            let r = io::wrap_decompress(r, in_path, compression)?;
            let fa = fasta::Reader::new(BufReader::new(r));
            let mut w = fasta::Writer::new(writer);
            for rec in fa.records() {
                let rec = rec?;
                processed += 1;
                if rng.random::<f64>() <= rate {
                    w.write(rec.id(), rec.desc(), rec.seq())?;
                    outputted += 1;
                }
                update_progress(pb);
            }
        }
        SeqFormat::Fastq => {
            let r = io::open_reader(in_path)?;
            let r = io::wrap_decompress(r, in_path, compression)?;
            let fq = fastq::Reader::new(BufReader::new(r));
            let mut w = fastq::Writer::new(writer);
            for rec in fq.records() {
                let rec = rec?;
                processed += 1;
                if rng.random::<f64>() <= rate {
                    w.write(rec.id(), rec.desc(), rec.seq(), rec.qual())?;
                    outputted += 1;
                }
                update_progress(pb);
            }
        }
        _ => bail!("sample currently supports FASTA/FASTQ input"),
    }
    finish_progress(pb);
    if !quiet {
        eprintln!("[INFO] {processed} sequences processed, {outputted} sequences outputted");
    }
    Ok(())
}

fn reservoir_from_iter<T>(
    items: impl Iterator<Item = T>,
    n: usize,
    rng: &mut ChaCha20Rng,
) -> Vec<T> {
    if n == 0 {
        return Vec::new();
    }
    let mut out = Vec::new();
    for (i, item) in items.enumerate() {
        if i < n {
            out.push(item);
        } else {
            let j = rng.random_range(0..=i);
            if j < n {
                out[j] = item;
            }
        }
    }
    out
}

fn reservoir_from_result_iter<T>(
    items: impl Iterator<Item = Result<T>>,
    n: usize,
    rng: &mut ChaCha20Rng,
    pb: &Option<ProgressBar>,
) -> Result<(Vec<T>, usize)> {
    if n == 0 {
        return Ok((Vec::new(), 0));
    }
    let mut out = Vec::new();
    let mut processed = 0usize;
    for (i, item) in items.enumerate() {
        let item = item?;
        processed += 1;
        if i < n {
            out.push(item);
        } else {
            let j = rng.random_range(0..=i);
            if j < n {
                out[j] = item;
            }
        }
        update_progress(pb);
    }
    Ok((out, processed))
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
