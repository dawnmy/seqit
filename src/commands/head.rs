use anyhow::{bail, Context, Result};
use bio::io::{fasta, fastq};
use std::io::BufReader;

use crate::{cli::HeadArgs, formats::SeqFormat, io, utils};

pub fn run(args: HeadArgs) -> Result<()> {
    let n = resolve_take(args.num, args.proportion)?;
    utils::validate_input_mode(
        "head",
        args.io.input.as_deref(),
        args.input1.as_deref(),
        args.input2.as_deref(),
    )?;
    if let (Some(in1), Some(in2)) = (args.input1.as_deref(), args.input2.as_deref()) {
        let out2 = args
            .output2
            .as_deref()
            .context("paired head requires --output2")?;
        if let Some(num) = args.num {
            return run_paired_streaming_num(&args, in1, in2, out2, num);
        }
        let r1 = io::read_records(Some(in1), SeqFormat::Fastq, &args.io.compression)?;
        let r2 = io::read_records(Some(in2), SeqFormat::Fastq, &args.io.compression)?;
        let keep = n(r1.len().min(r2.len()));
        io::write_records(
            &args.io.output,
            SeqFormat::Fastq,
            &args.io.compression,
            &r1[..keep],
        )?;
        io::write_records(out2, SeqFormat::Fastq, &args.io.compression, &r2[..keep])?;
        return Ok(());
    }

    let in_path = args.io.input.as_deref();
    if let Some(num) = args.num {
        return run_single_streaming_num(&args, in_path, num);
    }
    let (fmt, recs) = io::read_records_with_format(in_path, &args.io.format, &args.io.compression)?;
    let keep = n(recs.len());
    io::write_records(&args.io.output, fmt, &args.io.compression, &recs[..keep])
}

fn run_single_streaming_num(args: &HeadArgs, in_path: Option<&str>, n: usize) -> Result<()> {
    let (fmt, r) = io::open_seq_reader(in_path, &args.io.format, &args.io.compression)?;
    let w = io::open_writer(&args.io.output)?;
    let w = io::wrap_compress(w, &args.io.output, &args.io.compression)?;
    match fmt {
        SeqFormat::Fasta => {
            let fa = fasta::Reader::new(r);
            let mut out = fasta::Writer::new(w);
            for rec in fa.records().take(n) {
                let rec = rec?;
                out.write(rec.id(), rec.desc(), rec.seq())?;
            }
        }
        SeqFormat::Fastq => {
            let fq = fastq::Reader::new(r);
            let mut out = fastq::Writer::new(w);
            for rec in fq.records().take(n) {
                let rec = rec?;
                out.write(rec.id(), rec.desc(), rec.seq(), rec.qual())?;
            }
        }
        _ => bail!("head currently supports FASTA/FASTQ input"),
    }
    Ok(())
}

fn run_paired_streaming_num(
    args: &HeadArgs,
    in1: &str,
    in2: &str,
    out2: &str,
    n: usize,
) -> Result<()> {
    let r1 = io::open_reader(Some(in1))?;
    let r1 = io::wrap_decompress(r1, Some(in1), &args.io.compression)?;
    let r2 = io::open_reader(Some(in2))?;
    let r2 = io::wrap_decompress(r2, Some(in2), &args.io.compression)?;
    let fq1 = fastq::Reader::new(BufReader::new(r1));
    let fq2 = fastq::Reader::new(BufReader::new(r2));
    let mut it1 = fq1.records();
    let mut it2 = fq2.records();

    let w1 = io::open_writer(&args.io.output)?;
    let w1 = io::wrap_compress(w1, &args.io.output, &args.io.compression)?;
    let w2 = io::open_writer(out2)?;
    let w2 = io::wrap_compress(w2, out2, &args.io.compression)?;
    let mut w1 = fastq::Writer::new(w1);
    let mut w2 = fastq::Writer::new(w2);

    let mut kept = 0usize;
    let mut invalid = 0usize;
    while kept < n {
        let a = it1.next().transpose()?;
        let b = it2.next().transpose()?;
        match (a, b) {
            (Some(a), Some(b)) => {
                if pair_key(a.id()) != pair_key(b.id()) {
                    invalid += 1;
                    if !args.allow_unpaired {
                        bail!(
                            "paired inputs contain invalid/unpaired records (first mismatch: R1={} R2={}); rerun with --allow-unpaired to continue",
                            a.id(),
                            b.id()
                        );
                    }
                    continue;
                }
                w1.write(a.id(), a.desc(), a.seq(), a.qual())?;
                w2.write(b.id(), b.desc(), b.seq(), b.qual())?;
                kept += 1;
            }
            (Some(a), None) => {
                invalid += 1;
                if !args.allow_unpaired {
                    bail!(
                        "paired inputs contain unmatched record in R1 ({}); rerun with --allow-unpaired to continue",
                        a.id()
                    );
                }
            }
            (None, Some(b)) => {
                invalid += 1;
                if !args.allow_unpaired {
                    bail!(
                        "paired inputs contain unmatched record in R2 ({}); rerun with --allow-unpaired to continue",
                        b.id()
                    );
                }
            }
            (None, None) => break,
        }
    }
    if invalid > 0 && args.allow_unpaired {
        eprintln!("warning: skipped {invalid} invalid/unpaired paired records");
    }
    Ok(())
}

fn pair_key(id: &str) -> &str {
    id.strip_suffix("/1")
        .or_else(|| id.strip_suffix("/2"))
        .unwrap_or(id)
}

fn resolve_take(
    num: Option<usize>,
    proportion: Option<f64>,
) -> Result<Box<dyn Fn(usize) -> usize>> {
    match (num, proportion) {
        (Some(_), Some(_)) => bail!("choose only one of --num and --proportion"),
        (None, None) => bail!("one of --num or --proportion is required"),
        (Some(n), None) => Ok(Box::new(move |len| n.min(len))),
        (None, Some(p)) => {
            if !(0.0..=1.0).contains(&p) {
                bail!("--proportion must be in [0,1]");
            }
            Ok(Box::new(move |len| ((len as f64) * p).ceil() as usize))
        }
    }
}
