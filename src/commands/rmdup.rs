use std::collections::HashMap;

use anyhow::{bail, Context, Result};

use crate::{
    cli::{DupBy, RmdupArgs},
    formats::{SeqFormat, SeqRecord},
    io, pairs,
    utils::hash64,
};

pub fn run(args: RmdupArgs) -> Result<()> {
    if args.keep_first && args.keep_last {
        bail!("choose only one of --keep-first or --keep-last");
    }
    let paired_in1 = args.input1.as_deref().or(args.in1.as_deref());
    let paired_in2 = args.input2.as_deref().or(args.in2.as_deref());
    if let (Some(in1), Some(in2)) = (paired_in1, paired_in2) {
        let r1 = io::read_records(Some(in1), SeqFormat::Fastq, &args.io.compression)?;
        let r2 = io::read_records(Some(in2), SeqFormat::Fastq, &args.io.compression)?;
        let (r1, r2) = pairs::prepare_paired_records(r1, r2, args.allow_unpaired)?;
        let out2 = args
            .output2
            .as_deref()
            .context("paired rmdup requires --output2")?;
        let mut map: HashMap<u64, (SeqRecord, SeqRecord, usize)> = HashMap::new();
        for (a, b) in r1.into_iter().zip(r2.into_iter()) {
            let key = hash_pair(&a, &b, &args.by);
            let entry = map.entry(key).or_insert((a.clone(), b.clone(), 0));
            entry.2 += 1;
            if args.keep_last {
                *entry = (a, b, entry.2);
            }
        }
        let mut o1 = Vec::new();
        let mut o2 = Vec::new();
        for (_, (mut a, b, c)) in map {
            if args.mark_dup && c > 1 {
                a.id = format!("{};dup={}", a.id, c - 1);
            }
            o1.push(a);
            o2.push(b);
        }
        io::write_records(&args.io.output, SeqFormat::Fastq, &args.io.compression, &o1)?;
        io::write_records(out2, SeqFormat::Fastq, &args.io.compression, &o2)?;
        return Ok(());
    }

    let in_path = args.io.input.as_deref();
    let fmt = SeqFormat::from_arg(&args.io.format).unwrap_or(SeqFormat::detect(in_path)?);
    let recs = io::read_records(in_path, fmt, &args.io.compression)?;
    let mut map: HashMap<u64, (SeqRecord, usize)> = HashMap::new();
    for r in recs {
        let key = hash_single(&r, &args.by);
        let entry = map.entry(key).or_insert((r.clone(), 0));
        entry.1 += 1;
        if args.keep_last {
            *entry = (r, entry.1);
        }
    }
    let mut out = Vec::new();
    for (_, (mut r, c)) in map {
        if args.mark_dup && c > 1 {
            r.id = format!("{};dup={}", r.id, c - 1);
        }
        out.push(r);
    }
    io::write_records(&args.io.output, fmt, &args.io.compression, &out)
}

fn hash_single(r: &SeqRecord, by: &DupBy) -> u64 {
    match by {
        DupBy::Id => hash64(&r.id),
        DupBy::Seq => hash64(&r.seq),
        DupBy::Full => hash64(&(r.id.as_str(), &r.seq, r.qual.as_deref().unwrap_or_default())),
    }
}

fn hash_pair(a: &SeqRecord, b: &SeqRecord, by: &DupBy) -> u64 {
    match by {
        DupBy::Id => hash64(&(a.id.as_str(), b.id.as_str())),
        DupBy::Seq => hash64(&(&a.seq, &b.seq)),
        DupBy::Full => hash64(&(
            a.id.as_str(),
            b.id.as_str(),
            &a.seq,
            &b.seq,
            a.qual.as_deref().unwrap_or_default(),
            b.qual.as_deref().unwrap_or_default(),
        )),
    }
}
