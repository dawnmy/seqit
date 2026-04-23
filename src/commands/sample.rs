use anyhow::{bail, Context, Result};
use bio::io::{fasta, fastq};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;
use std::io::{BufReader, BufWriter};

use crate::{
    cli::SampleArgs,
    formats::{SeqFormat, SeqRecord},
    io,
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
        let out2 = args
            .output2
            .as_deref()
            .context("paired sample requires --output2")?;
        if let Some(rate) = args.rate {
            log_info(&args, "sample by rate");
            let outputted = sample_paired_rate_streaming(
                in1,
                in2,
                &args.io.output,
                out2,
                &args.io.compression,
                rate,
                &mut rng,
                args.allow_unpaired,
            )?;
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
        io::write_records(&args.io.output, SeqFormat::Fastq, &args.io.compression, &o1)?;
        io::write_records(out2, SeqFormat::Fastq, &args.io.compression, &o2)?;
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
        return sample_rate_streaming(
            reader,
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
            let fa = fasta::Reader::new(reader);
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
            let fq = fastq::Reader::new(reader);
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

fn sample_paired_rate_streaming(
    in1: &str,
    in2: &str,
    out1: &str,
    out2: &str,
    compression: &crate::cli::CompressionArg,
    rate: f64,
    rng: &mut ChaCha20Rng,
    allow_unpaired: bool,
) -> Result<usize> {
    let r1 = io::open_reader(Some(in1))?;
    let r1 = io::wrap_decompress(r1, Some(in1), compression)?;
    let r2 = io::open_reader(Some(in2))?;
    let r2 = io::wrap_decompress(r2, Some(in2), compression)?;
    let fq1 = fastq::Reader::new(BufReader::new(r1));
    let fq2 = fastq::Reader::new(BufReader::new(r2));
    let mut it1 = fq1.records();
    let mut it2 = fq2.records();

    let w1 = io::open_writer(out1)?;
    let w1 = io::wrap_compress(w1, out1, compression)?;
    let w2 = io::open_writer(out2)?;
    let w2 = io::wrap_compress(w2, out2, compression)?;
    let mut w1 = fastq::Writer::new(BufWriter::new(w1));
    let mut w2 = fastq::Writer::new(BufWriter::new(w2));

    let mut outputted = 0usize;
    let mut invalid_preview = Vec::new();
    let mut invalid_count = 0usize;
    loop {
        let a = it1.next().transpose()?;
        let b = it2.next().transpose()?;
        match (a, b) {
            (Some(a), Some(b)) => {
                if pair_key(a.id()) != pair_key(b.id()) {
                    invalid_count += 1;
                    if invalid_preview.len() < 10 {
                        invalid_preview.push(format!("R1={} R2={}", a.id(), b.id()));
                    }
                    continue;
                }
                if rng.random::<f64>() <= rate {
                    w1.write(a.id(), a.desc(), a.seq(), a.qual())?;
                    w2.write(b.id(), b.desc(), b.seq(), b.qual())?;
                    outputted += 1;
                }
            }
            (Some(a), None) => {
                invalid_count += 1;
                if invalid_preview.len() < 10 {
                    invalid_preview.push(format!("unmatched in R1: {}", a.id()));
                }
            }
            (None, Some(b)) => {
                invalid_count += 1;
                if invalid_preview.len() < 10 {
                    invalid_preview.push(format!("unmatched in R2: {}", b.id()));
                }
            }
            (None, None) => break,
        }
    }
    report_invalid_pairs(invalid_count, invalid_preview, allow_unpaired)?;
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
    let r1 = io::open_reader(Some(in1))?;
    let r1 = io::wrap_decompress(r1, Some(in1), compression)?;
    let r2 = io::open_reader(Some(in2))?;
    let r2 = io::wrap_decompress(r2, Some(in2), compression)?;
    let fq1 = fastq::Reader::new(BufReader::new(r1));
    let fq2 = fastq::Reader::new(BufReader::new(r2));
    let mut it1 = fq1.records();
    let mut it2 = fq2.records();

    let mut output = Vec::with_capacity(n);
    let mut valid_seen = 0usize;
    let mut invalid_preview = Vec::new();
    let mut invalid_count = 0usize;
    loop {
        let a = it1.next().transpose()?;
        let b = it2.next().transpose()?;
        match (a, b) {
            (Some(a), Some(b)) => {
                if pair_key(a.id()) != pair_key(b.id()) {
                    invalid_count += 1;
                    if invalid_preview.len() < 10 {
                        invalid_preview.push(format!("R1={} R2={}", a.id(), b.id()));
                    }
                    continue;
                }
                valid_seen += 1;
                let pair = (
                    SeqRecord {
                        id: a.id().to_string(),
                        desc: a.desc().map(|d| d.to_string()),
                        seq: a.seq().to_vec(),
                        qual: Some(a.qual().to_vec()),
                    },
                    SeqRecord {
                        id: b.id().to_string(),
                        desc: b.desc().map(|d| d.to_string()),
                        seq: b.seq().to_vec(),
                        qual: Some(b.qual().to_vec()),
                    },
                );
                if output.len() < n {
                    output.push(pair);
                } else if n > 0 {
                    let j = rng.random_range(0..valid_seen);
                    if j < n {
                        output[j] = pair;
                    }
                }
            }
            (Some(a), None) => {
                invalid_count += 1;
                if invalid_preview.len() < 10 {
                    invalid_preview.push(format!("unmatched in R1: {}", a.id()));
                }
            }
            (None, Some(b)) => {
                invalid_count += 1;
                if invalid_preview.len() < 10 {
                    invalid_preview.push(format!("unmatched in R2: {}", b.id()));
                }
            }
            (None, None) => break,
        }
    }
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
    reader: BufReader<Box<dyn std::io::Read>>,
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
            let fa = fasta::Reader::new(reader);
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
            let fq = fastq::Reader::new(reader);
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

fn reservoir_from_result_iter<T>(
    items: impl Iterator<Item = Result<T>>,
    n: usize,
    rng: &mut ChaCha20Rng,
    pb: &Option<ProgressBar>,
) -> Result<(Vec<T>, usize)> {
    if n == 0 {
        return Ok((Vec::new(), 0));
    }
    let mut items = items;
    let mut out = Vec::with_capacity(n);
    let mut processed = 0usize;

    for _ in 0..n {
        let Some(item) = items.next() else {
            return Ok((out, processed));
        };
        let item = item?;
        processed += 1;
        out.push(item);
        update_progress(pb);
    }

    // Vitter's Algorithm L (skip-based reservoir sampling).
    // This greatly reduces RNG calls on huge inputs where n << total records.
    let mut w = (sample_open01(rng).ln() / n as f64).exp();
    loop {
        let gap = (sample_open01(rng).ln() / (1.0 - w).ln()).floor() as usize;
        for _ in 0..gap {
            let Some(item) = items.next() else {
                return Ok((out, processed));
            };
            item?;
            processed += 1;
            update_progress(pb);
        }

        let Some(item) = items.next() else {
            return Ok((out, processed));
        };
        let item = item?;
        processed += 1;
        update_progress(pb);

        let j = rng.random_range(0..n);
        out[j] = item;

        w *= (sample_open01(rng).ln() / n as f64).exp();
    }
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
