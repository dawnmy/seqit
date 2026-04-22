use anyhow::{bail, Context, Result};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;

use crate::{cli::SpikeArgs, formats::SeqFormat, io, pairs};

pub fn run(args: SpikeArgs) -> Result<()> {
    let paired_in1 = args.input1.as_deref().or(args.in1.as_deref());
    let paired_in2 = args.input2.as_deref().or(args.in2.as_deref());
    if let (Some(in1), Some(in2), Some(add2)) = (paired_in1, paired_in2, args.add2.as_deref()) {
        let t1 = io::read_records(Some(in1), SeqFormat::Fastq, &args.compression)?;
        let t2 = io::read_records(Some(in2), SeqFormat::Fastq, &args.compression)?;
        let a1 = io::read_records(Some(&args.add), SeqFormat::Fastq, &args.compression)?;
        let a2 = io::read_records(Some(add2), SeqFormat::Fastq, &args.compression)?;
        let (t1, t2) = pairs::prepare_paired_records(t1, t2, args.allow_unpaired)?;
        let (a1, a2) = pairs::prepare_paired_records(a1, a2, args.allow_unpaired)?;
        let out2 = args
            .output2
            .as_deref()
            .context("paired spike requires --output2")?;
        let (o1, o2) = spike_pairs(t1, t2, a1, a2, args.seed);
        io::write_records(&args.output, SeqFormat::Fastq, &args.compression, &o1)?;
        io::write_records(out2, SeqFormat::Fastq, &args.compression, &o2)?;
        return Ok(());
    }

    if paired_in1.is_some() || paired_in2.is_some() || args.add2.is_some() {
        bail!("paired spike requires -i/-I (or --in1/--in2) and -a/-A together");
    }
    let input = args.input.as_deref();
    let fmt = SeqFormat::from_arg(&args.format).unwrap_or(SeqFormat::detect(input)?);
    let target = io::read_records(input, fmt, &args.compression)?;
    let add = io::read_records(Some(&args.add), fmt, &args.compression)?;
    let out = spike_single(target, add, args.seed);
    io::write_records(&args.output, fmt, &args.compression, &out)
}

fn sample_slots(n_slots: usize, m: usize, seed: u64) -> Vec<usize> {
    let mut rng = ChaCha20Rng::seed_from_u64(seed);
    (0..m).map(|_| rng.random_range(0..n_slots)).collect()
}

fn spike_single(
    target: Vec<crate::formats::SeqRecord>,
    insert: Vec<crate::formats::SeqRecord>,
    seed: u64,
) -> Vec<crate::formats::SeqRecord> {
    let mut buckets = vec![Vec::new(); target.len() + 1];
    let slots = sample_slots(target.len() + 1, insert.len(), seed);
    for (ins, slot) in insert.into_iter().zip(slots) {
        buckets[slot].push(ins);
    }
    let mut out = Vec::new();
    out.append(&mut buckets[0]);
    for (i, t) in target.into_iter().enumerate() {
        out.push(t);
        out.append(&mut buckets[i + 1]);
    }
    out
}

fn spike_pairs(
    t1: Vec<crate::formats::SeqRecord>,
    t2: Vec<crate::formats::SeqRecord>,
    a1: Vec<crate::formats::SeqRecord>,
    a2: Vec<crate::formats::SeqRecord>,
    seed: u64,
) -> (
    Vec<crate::formats::SeqRecord>,
    Vec<crate::formats::SeqRecord>,
) {
    let mut b1 = vec![Vec::new(); t1.len() + 1];
    let mut b2 = vec![Vec::new(); t1.len() + 1];
    let slots = sample_slots(t1.len() + 1, a1.len(), seed);
    for ((x1, x2), slot) in a1.into_iter().zip(a2.into_iter()).zip(slots) {
        b1[slot].push(x1);
        b2[slot].push(x2);
    }
    let mut o1 = Vec::new();
    let mut o2 = Vec::new();
    o1.append(&mut b1[0]);
    o2.append(&mut b2[0]);
    for (i, (x1, x2)) in t1.into_iter().zip(t2.into_iter()).enumerate() {
        o1.push(x1);
        o2.push(x2);
        o1.append(&mut b1[i + 1]);
        o2.append(&mut b2[i + 1]);
    }
    (o1, o2)
}
