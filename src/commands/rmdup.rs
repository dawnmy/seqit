use ahash::{AHashMap, AHashSet};
use anyhow::{bail, Context, Result};
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
    io::for_each_record_from_reader(br, fmt, |rec| {
        if seen.insert(single_key_ref(&rec, &args.by)?) {
            io::write_record_ref(&mut out, fmt, &rec)?;
        }
        Ok(())
    })?;
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
    let w1 = io::open_writer(&args.io.output)?;
    let w1 = io::wrap_compress_for_streams(w1, &args.io.output, &args.io.compression, 2)?;
    let w2 = io::open_writer(out2)?;
    let w2 = io::wrap_compress_for_streams(w2, out2, &args.io.compression, 2)?;
    let mut w1 = io::buffered_writer(w1);
    let mut w2 = io::buffered_writer(w2);

    let mut seen = AHashSet::new();
    let mut invalid_preview = Vec::new();
    let mut invalid_count = 0usize;
    io::for_each_fastq_pair(in1, in2, &args.io.compression, |pair| {
        match pair {
            io::FastqPair::Both(a, b) => {
                let a_id = a.id_str()?;
                let b_id = b.id_str()?;
                if pair_id_key(a_id) != pair_id_key(b_id) {
                    invalid_count += 1;
                    if invalid_preview.len() < 10 {
                        invalid_preview.push(format!("R1={a_id} R2={b_id}"));
                    }
                    return Ok(true);
                }
                if seen.insert(pair_key_ref(&a, &b, &args.by)?) {
                    io::write_record_ref(&mut w1, SeqFormat::Fastq, &a)?;
                    io::write_record_ref(&mut w2, SeqFormat::Fastq, &b)?;
                }
            }
            io::FastqPair::Left(a) => {
                record_unmatched(&mut invalid_count, &mut invalid_preview, "R1", a.id_str()?)
            }
            io::FastqPair::Right(b) => {
                record_unmatched(&mut invalid_count, &mut invalid_preview, "R2", b.id_str()?)
            }
        }
        Ok(true)
    })?;
    report_invalid_pairs(invalid_count, invalid_preview, args.allow_unpaired)
}

fn run_paired_map(args: &RmdupArgs, in1: &str, in2: &str, out2: &str) -> Result<()> {
    let mut map: AHashMap<PairKey, (SeqRecord, SeqRecord, usize)> = AHashMap::new();
    let mut invalid_preview = Vec::new();
    let mut invalid_count = 0usize;
    io::for_each_fastq_pair(in1, in2, &args.io.compression, |pair| {
        match pair {
            io::FastqPair::Both(a, b) => {
                let a_id = a.id_str()?;
                let b_id = b.id_str()?;
                if pair_id_key(a_id) != pair_id_key(b_id) {
                    invalid_count += 1;
                    if invalid_preview.len() < 10 {
                        invalid_preview.push(format!("R1={a_id} R2={b_id}"));
                    }
                    return Ok(true);
                }
                let key = pair_key_ref(&a, &b, &args.by)?;
                match map.entry(key) {
                    Entry::Vacant(entry) => {
                        let left = a.to_owned_record()?;
                        let right = b.to_owned_record()?;
                        entry.insert((left, right, 1));
                    }
                    Entry::Occupied(mut entry) => {
                        let value = entry.get_mut();
                        value.2 += 1;
                        if args.keep_last {
                            value.0 = a.to_owned_record()?;
                            value.1 = b.to_owned_record()?;
                        }
                    }
                }
            }
            io::FastqPair::Left(a) => {
                record_unmatched(&mut invalid_count, &mut invalid_preview, "R1", a.id_str()?)
            }
            io::FastqPair::Right(b) => {
                record_unmatched(&mut invalid_count, &mut invalid_preview, "R2", b.id_str()?)
            }
        }
        Ok(true)
    })?;
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

fn single_key_ref(r: &io::SeqRecordRef<'_>, by: &DupBy) -> Result<SingleKey> {
    Ok(match by {
        DupBy::Id => SingleKey::Id(r.id_str()?.to_string()),
        DupBy::Seq => SingleKey::Seq(r.seq.to_vec()),
        DupBy::Full => SingleKey::Full(
            r.id_str()?.to_string(),
            r.seq.to_vec(),
            r.qual.map(<[u8]>::to_vec),
        ),
    })
}

fn pair_key_ref(a: &io::SeqRecordRef<'_>, b: &io::SeqRecordRef<'_>, by: &DupBy) -> Result<PairKey> {
    Ok(match by {
        DupBy::Id => PairKey::Id(a.id_str()?.to_string(), b.id_str()?.to_string()),
        DupBy::Seq => PairKey::Seq(a.seq.to_vec(), b.seq.to_vec()),
        DupBy::Full => PairKey::Full(
            a.id_str()?.to_string(),
            b.id_str()?.to_string(),
            a.seq.to_vec(),
            b.seq.to_vec(),
            a.qual.map(<[u8]>::to_vec),
            b.qual.map(<[u8]>::to_vec),
        ),
    })
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
