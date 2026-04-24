use ahash::{AHashMap, AHashSet};
use anyhow::{bail, Context, Result};
use bio::io::{fasta, fastq};
use std::collections::hash_map::Entry;

use crate::{
    cli::{DupBy, RmdupArgs},
    formats::{SeqFormat, SeqRecord},
    io, utils,
};

pub fn run(args: RmdupArgs) -> Result<()> {
    if args.keep_first && args.keep_last {
        bail!("choose only one of --keep-first or --keep-last");
    }
    let paired_in1 = args.input1.as_deref().or(args.in1.as_deref());
    let paired_in2 = args.input2.as_deref().or(args.in2.as_deref());
    utils::validate_input_mode("rmdup", args.io.input.as_deref(), paired_in1, paired_in2)?;
    if let (Some(in1), Some(in2)) = (paired_in1, paired_in2) {
        let out2 = args
            .output2
            .as_deref()
            .context("paired rmdup requires --output2")?;
        if !args.keep_last && !args.mark_dup && !args.count_dup {
            return run_paired_keep_first_streaming(&args, in1, in2, out2);
        }
        return run_paired_map(&args, in1, in2, out2);
    }

    let in_path = args.io.input.as_deref();
    if !args.keep_last && !args.mark_dup && !args.count_dup {
        return run_single_keep_first_streaming(&args, in_path);
    }
    run_single_map(&args, in_path)
}

fn run_single_keep_first_streaming(args: &RmdupArgs, in_path: Option<&str>) -> Result<()> {
    let (fmt, br) = io::open_seq_reader(in_path, &args.io.format, &args.io.compression)?;
    let w = io::open_writer(&args.io.output)?;
    let w = io::wrap_compress(w, &args.io.output, &args.io.compression)?;
    let mut out = io::buffered_writer(w);
    let mut seen = AHashSet::new();
    match fmt {
        SeqFormat::Fasta => {
            let fa = fasta::Reader::new(br);
            for rec in fa.records() {
                let rec = rec?;
                let seq_rec = SeqRecord {
                    id: rec.id().to_string(),
                    desc: rec.desc().map(|d| d.to_string()),
                    seq: rec.seq().to_vec(),
                    qual: None,
                };
                if seen.insert(single_key(&seq_rec, &args.by)) {
                    io::write_record(&mut out, fmt, &seq_rec)?;
                }
            }
        }
        SeqFormat::Fastq => {
            let fq = fastq::Reader::new(br);
            for rec in fq.records() {
                let rec = rec?;
                let seq_rec = SeqRecord {
                    id: rec.id().to_string(),
                    desc: rec.desc().map(|d| d.to_string()),
                    seq: rec.seq().to_vec(),
                    qual: Some(rec.qual().to_vec()),
                };
                if seen.insert(single_key(&seq_rec, &args.by)) {
                    io::write_record(&mut out, fmt, &seq_rec)?;
                }
            }
        }
        _ => bail!("rmdup currently supports FASTA/FASTQ input"),
    }
    Ok(())
}

fn run_single_map(args: &RmdupArgs, in_path: Option<&str>) -> Result<()> {
    let mut map: AHashMap<SingleKey, (SeqRecord, usize)> = AHashMap::new();
    let fmt = io::for_each_record(in_path, &args.io.format, &args.io.compression, |r| {
        let key = single_key(&r, &args.by);
        match map.entry(key) {
            Entry::Vacant(entry) => {
                entry.insert((r, 1));
            }
            Entry::Occupied(mut entry) => {
                let value = entry.get_mut();
                value.1 += 1;
                if args.keep_last {
                    value.0 = r;
                }
            }
        }
        Ok(())
    })?;
    let mut out = Vec::new();
    for (_, (mut r, c)) in map {
        annotate_id(&mut r.id, c, args);
        out.push(r);
    }
    io::write_records(&args.io.output, fmt, &args.io.compression, &out)
}

fn run_paired_keep_first_streaming(
    args: &RmdupArgs,
    in1: &str,
    in2: &str,
    out2: &str,
) -> Result<()> {
    let r1 = io::open_reader(Some(in1))?;
    let r1 = io::wrap_decompress(r1, Some(in1), &args.io.compression)?;
    let r2 = io::open_reader(Some(in2))?;
    let r2 = io::wrap_decompress(r2, Some(in2), &args.io.compression)?;
    let fq1 = fastq::Reader::new(io::buffered_reader(r1));
    let fq2 = fastq::Reader::new(io::buffered_reader(r2));
    let mut it1 = fq1.records();
    let mut it2 = fq2.records();

    let w1 = io::open_writer(&args.io.output)?;
    let w1 = io::wrap_compress_for_streams(w1, &args.io.output, &args.io.compression, 2)?;
    let w2 = io::open_writer(out2)?;
    let w2 = io::wrap_compress_for_streams(w2, out2, &args.io.compression, 2)?;
    let mut w1 = fastq::Writer::new(io::buffered_writer(w1));
    let mut w2 = fastq::Writer::new(io::buffered_writer(w2));

    let mut seen = AHashSet::new();
    let mut invalid_preview = Vec::new();
    let mut invalid_count = 0usize;
    loop {
        let a = it1.next().transpose()?;
        let b = it2.next().transpose()?;
        match (a, b) {
            (Some(a), Some(b)) => {
                if pair_id_key(a.id()) != pair_id_key(b.id()) {
                    invalid_count += 1;
                    if invalid_preview.len() < 10 {
                        invalid_preview.push(format!("R1={} R2={}", a.id(), b.id()));
                    }
                    continue;
                }
                if seen.insert(pair_key_from_fastq(&a, &b, &args.by)) {
                    w1.write(a.id(), a.desc(), a.seq(), a.qual())?;
                    w2.write(b.id(), b.desc(), b.seq(), b.qual())?;
                }
            }
            (Some(a), None) => {
                record_unmatched(&mut invalid_count, &mut invalid_preview, "R1", a.id())
            }
            (None, Some(b)) => {
                record_unmatched(&mut invalid_count, &mut invalid_preview, "R2", b.id())
            }
            (None, None) => break,
        }
    }
    report_invalid_pairs(invalid_count, invalid_preview, args.allow_unpaired)
}

fn run_paired_map(args: &RmdupArgs, in1: &str, in2: &str, out2: &str) -> Result<()> {
    let r1 = io::open_reader(Some(in1))?;
    let r1 = io::wrap_decompress(r1, Some(in1), &args.io.compression)?;
    let r2 = io::open_reader(Some(in2))?;
    let r2 = io::wrap_decompress(r2, Some(in2), &args.io.compression)?;
    let fq1 = fastq::Reader::new(io::buffered_reader(r1));
    let fq2 = fastq::Reader::new(io::buffered_reader(r2));
    let mut it1 = fq1.records();
    let mut it2 = fq2.records();

    let mut map: AHashMap<PairKey, (SeqRecord, SeqRecord, usize)> = AHashMap::new();
    let mut invalid_preview = Vec::new();
    let mut invalid_count = 0usize;
    loop {
        let a = it1.next().transpose()?;
        let b = it2.next().transpose()?;
        match (a, b) {
            (Some(a), Some(b)) => {
                if pair_id_key(a.id()) != pair_id_key(b.id()) {
                    invalid_count += 1;
                    if invalid_preview.len() < 10 {
                        invalid_preview.push(format!("R1={} R2={}", a.id(), b.id()));
                    }
                    continue;
                }
                let left = SeqRecord {
                    id: a.id().to_string(),
                    desc: a.desc().map(|d| d.to_string()),
                    seq: a.seq().to_vec(),
                    qual: Some(a.qual().to_vec()),
                };
                let right = SeqRecord {
                    id: b.id().to_string(),
                    desc: b.desc().map(|d| d.to_string()),
                    seq: b.seq().to_vec(),
                    qual: Some(b.qual().to_vec()),
                };
                let key = pair_key(&left, &right, &args.by);
                match map.entry(key) {
                    Entry::Vacant(entry) => {
                        entry.insert((left, right, 1));
                    }
                    Entry::Occupied(mut entry) => {
                        let value = entry.get_mut();
                        value.2 += 1;
                        if args.keep_last {
                            value.0 = left;
                            value.1 = right;
                        }
                    }
                }
            }
            (Some(a), None) => {
                record_unmatched(&mut invalid_count, &mut invalid_preview, "R1", a.id())
            }
            (None, Some(b)) => {
                record_unmatched(&mut invalid_count, &mut invalid_preview, "R2", b.id())
            }
            (None, None) => break,
        }
    }
    report_invalid_pairs(invalid_count, invalid_preview, args.allow_unpaired)?;

    let mut o1 = Vec::with_capacity(map.len());
    let mut o2 = Vec::with_capacity(map.len());
    for (_, (mut a, b, c)) in map {
        annotate_id(&mut a.id, c, args);
        o1.push(a);
        o2.push(b);
    }
    io::write_record_pair_parallel(
        &args.io.output,
        out2,
        SeqFormat::Fastq,
        &args.io.compression,
        &o1,
        &o2,
    )
}

#[derive(Hash, PartialEq, Eq)]
enum SingleKey {
    Id(String),
    Seq(Vec<u8>),
    Full(String, Vec<u8>, Option<Vec<u8>>),
}

#[derive(Hash, PartialEq, Eq)]
enum PairKey {
    Id(String, String),
    Seq(Vec<u8>, Vec<u8>),
    Full(
        String,
        String,
        Vec<u8>,
        Vec<u8>,
        Option<Vec<u8>>,
        Option<Vec<u8>>,
    ),
}

fn single_key(r: &SeqRecord, by: &DupBy) -> SingleKey {
    match by {
        DupBy::Id => SingleKey::Id(r.id.clone()),
        DupBy::Seq => SingleKey::Seq(r.seq.clone()),
        DupBy::Full => SingleKey::Full(r.id.clone(), r.seq.clone(), r.qual.clone()),
    }
}

fn pair_key(a: &SeqRecord, b: &SeqRecord, by: &DupBy) -> PairKey {
    match by {
        DupBy::Id => PairKey::Id(a.id.clone(), b.id.clone()),
        DupBy::Seq => PairKey::Seq(a.seq.clone(), b.seq.clone()),
        DupBy::Full => PairKey::Full(
            a.id.clone(),
            b.id.clone(),
            a.seq.clone(),
            b.seq.clone(),
            a.qual.clone(),
            b.qual.clone(),
        ),
    }
}

fn pair_key_from_fastq(a: &fastq::Record, b: &fastq::Record, by: &DupBy) -> PairKey {
    match by {
        DupBy::Id => PairKey::Id(a.id().to_string(), b.id().to_string()),
        DupBy::Seq => PairKey::Seq(a.seq().to_vec(), b.seq().to_vec()),
        DupBy::Full => PairKey::Full(
            a.id().to_string(),
            b.id().to_string(),
            a.seq().to_vec(),
            b.seq().to_vec(),
            Some(a.qual().to_vec()),
            Some(b.qual().to_vec()),
        ),
    }
}

fn annotate_id(id: &mut String, count: usize, args: &RmdupArgs) {
    if args.count_dup {
        id.push_str(&format!(";count={count}"));
    }
    if args.mark_dup && count > 1 {
        id.push_str(&format!(";dup={}", count - 1));
    }
}

fn pair_id_key(id: &str) -> &str {
    id.strip_suffix("/1")
        .or_else(|| id.strip_suffix("/2"))
        .unwrap_or(id)
}

fn record_unmatched(count: &mut usize, preview: &mut Vec<String>, mate: &str, id: &str) {
    *count += 1;
    if preview.len() < 10 {
        preview.push(format!("unmatched in {mate}: {id}"));
    }
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
