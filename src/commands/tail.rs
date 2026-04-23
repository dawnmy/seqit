use anyhow::{bail, Context, Result};
use bio::io::{fasta, fastq};
use std::collections::VecDeque;
use std::io::{BufReader, BufWriter};

use crate::{cli::TailArgs, formats::SeqFormat, io, utils};

pub fn run(args: TailArgs) -> Result<()> {
    let _ = resolve_take(args.num, args.proportion)?;
    utils::validate_input_mode(
        "tail",
        args.io.input.as_deref(),
        args.input1.as_deref(),
        args.input2.as_deref(),
    )?;
    if let (Some(in1), Some(in2)) = (args.input1.as_deref(), args.input2.as_deref()) {
        let out2 = args
            .output2
            .as_deref()
            .context("paired tail requires --output2")?;
        return run_paired_tail_streaming(&args, in1, in2, out2);
    }

    let in_path = args.io.input.as_deref();
    run_single_tail_streaming(&args, in_path)
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

fn run_single_tail_streaming(args: &TailArgs, in_path: Option<&str>) -> Result<()> {
    let keep = resolve_keep_count_single(args, in_path)?;
    let (fmt, br) = io::open_seq_reader(in_path, &args.io.format, &args.io.compression)?;
    let w = io::open_writer(&args.io.output)?;
    let w = io::wrap_compress(w, &args.io.output, &args.io.compression)?;
    let mut bw = BufWriter::new(w);

    match fmt {
        SeqFormat::Fasta => {
            let fa = fasta::Reader::new(br);
            let mut ring: VecDeque<(String, Option<String>, Vec<u8>)> =
                VecDeque::with_capacity(keep);
            for rec in fa.records() {
                let rec = rec?;
                if keep == 0 {
                    continue;
                }
                if ring.len() == keep {
                    ring.pop_front();
                }
                ring.push_back((
                    rec.id().to_string(),
                    rec.desc().map(|d| d.to_string()),
                    rec.seq().to_vec(),
                ));
            }
            let mut out = fasta::Writer::new(&mut bw);
            for (id, desc, seq) in ring {
                out.write(&id, desc.as_deref(), &seq)?;
            }
        }
        SeqFormat::Fastq => {
            let fq = fastq::Reader::new(br);
            let mut ring: VecDeque<(String, Option<String>, Vec<u8>, Vec<u8>)> =
                VecDeque::with_capacity(keep);
            for rec in fq.records() {
                let rec = rec?;
                if keep == 0 {
                    continue;
                }
                if ring.len() == keep {
                    ring.pop_front();
                }
                ring.push_back((
                    rec.id().to_string(),
                    rec.desc().map(|d| d.to_string()),
                    rec.seq().to_vec(),
                    rec.qual().to_vec(),
                ));
            }
            let mut out = fastq::Writer::new(&mut bw);
            for (id, desc, seq, qual) in ring {
                out.write(&id, desc.as_deref(), &seq, &qual)?;
            }
        }
        _ => bail!("tail currently supports FASTA/FASTQ input"),
    }
    Ok(())
}

fn resolve_keep_count_single(args: &TailArgs, in_path: Option<&str>) -> Result<usize> {
    match (args.num, args.proportion) {
        (Some(n), None) => Ok(n),
        (None, Some(p)) => {
            if in_path == Some("-") {
                // stdin cannot be rewound; fallback to materialize
                let (_fmt, recs) =
                    io::read_records_with_format(in_path, &args.io.format, &args.io.compression)?;
                return Ok(((recs.len() as f64) * p).ceil() as usize);
            }
            let (fmt, br) = io::open_seq_reader(in_path, &args.io.format, &args.io.compression)?;
            let total = count_records(fmt, br)?;
            Ok(((total as f64) * p).ceil() as usize)
        }
        _ => unreachable!(),
    }
}

fn count_records(fmt: SeqFormat, br: BufReader<Box<dyn std::io::Read>>) -> Result<usize> {
    let mut total = 0usize;
    match fmt {
        SeqFormat::Fasta => {
            let fa = fasta::Reader::new(br);
            for rec in fa.records() {
                rec?;
                total += 1;
            }
        }
        SeqFormat::Fastq => {
            let fq = fastq::Reader::new(br);
            for rec in fq.records() {
                rec?;
                total += 1;
            }
        }
        _ => bail!("tail currently supports FASTA/FASTQ input"),
    }
    Ok(total)
}

fn run_paired_tail_streaming(args: &TailArgs, in1: &str, in2: &str, out2: &str) -> Result<()> {
    let keep = resolve_keep_count_paired(args, in1, in2)?;
    let r1 = io::open_reader(Some(in1))?;
    let r1 = io::wrap_decompress(r1, Some(in1), &args.io.compression)?;
    let r2 = io::open_reader(Some(in2))?;
    let r2 = io::wrap_decompress(r2, Some(in2), &args.io.compression)?;
    let fq1 = fastq::Reader::new(BufReader::new(r1));
    let fq2 = fastq::Reader::new(BufReader::new(r2));
    let mut it1 = fq1.records();
    let mut it2 = fq2.records();

    let mut ring: VecDeque<(
        (String, Option<String>, Vec<u8>, Vec<u8>),
        (String, Option<String>, Vec<u8>, Vec<u8>),
    )> = VecDeque::with_capacity(keep);
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
                if keep == 0 {
                    continue;
                }
                if ring.len() == keep {
                    ring.pop_front();
                }
                ring.push_back((
                    (
                        a.id().to_string(),
                        a.desc().map(|d| d.to_string()),
                        a.seq().to_vec(),
                        a.qual().to_vec(),
                    ),
                    (
                        b.id().to_string(),
                        b.desc().map(|d| d.to_string()),
                        b.seq().to_vec(),
                        b.qual().to_vec(),
                    ),
                ));
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
    report_invalid_pairs(invalid_count, invalid_preview, args.allow_unpaired)?;

    let w1 = io::open_writer(&args.io.output)?;
    let w1 = io::wrap_compress(w1, &args.io.output, &args.io.compression)?;
    let w2 = io::open_writer(out2)?;
    let w2 = io::wrap_compress(w2, out2, &args.io.compression)?;
    let mut w1 = fastq::Writer::new(BufWriter::new(w1));
    let mut w2 = fastq::Writer::new(BufWriter::new(w2));
    for ((id1, desc1, seq1, q1), (id2, desc2, seq2, q2)) in ring {
        w1.write(&id1, desc1.as_deref(), &seq1, &q1)?;
        w2.write(&id2, desc2.as_deref(), &seq2, &q2)?;
    }
    Ok(())
}

fn resolve_keep_count_paired(args: &TailArgs, in1: &str, in2: &str) -> Result<usize> {
    match (args.num, args.proportion) {
        (Some(n), None) => Ok(n),
        (None, Some(p)) => {
            let r1 = io::open_reader(Some(in1))?;
            let r1 = io::wrap_decompress(r1, Some(in1), &args.io.compression)?;
            let r2 = io::open_reader(Some(in2))?;
            let r2 = io::wrap_decompress(r2, Some(in2), &args.io.compression)?;
            let fq1 = fastq::Reader::new(BufReader::new(r1));
            let fq2 = fastq::Reader::new(BufReader::new(r2));
            let mut it1 = fq1.records();
            let mut it2 = fq2.records();
            let mut valid = 0usize;
            loop {
                let a = it1.next().transpose()?;
                let b = it2.next().transpose()?;
                match (a, b) {
                    (Some(a), Some(b)) => {
                        if pair_key(a.id()) == pair_key(b.id()) {
                            valid += 1;
                        }
                    }
                    (None, None) => break,
                    _ => {}
                }
            }
            Ok(((valid as f64) * p).ceil() as usize)
        }
        _ => unreachable!(),
    }
}

fn pair_key(id: &str) -> &str {
    id.strip_suffix("/1")
        .or_else(|| id.strip_suffix("/2"))
        .unwrap_or(id)
}

fn report_invalid_pairs(
    invalid_count: usize,
    invalid_preview: Vec<String>,
    allow_unpaired: bool,
) -> Result<()> {
    if invalid_count == 0 {
        return Ok(());
    }
    let preview = invalid_preview.join("; ");
    if !allow_unpaired {
        bail!(
            "paired inputs contain {invalid_count} invalid/unpaired records; examples: {preview}. rerun with --allow-unpaired to continue"
        );
    }
    eprintln!(
        "warning: paired inputs contain {invalid_count} invalid/unpaired records; examples: {preview}. unpaired records will be skipped"
    );
    Ok(())
}
