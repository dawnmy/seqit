use anyhow::{bail, Context, Result};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;

use crate::{cli::SpikeArgs, formats::SeqFormat, io, pairs, utils};

pub fn run(args: SpikeArgs) -> Result<()> {
    let paired_in2 = args.input2.as_deref().or(args.in2.as_deref());
    let paired_in1 = args.input1.as_deref().or(args.in1.as_deref());
    let single_input = args.input.as_deref().or(args.input_legacy.as_deref());
    utils::validate_input_mode("spike", single_input, paired_in1, paired_in2)?;

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
        bail!(
            "paired spike requires R1 + R2 target inputs and R1 + R2 spike-in inputs together (use -i/-I with -a/-A, or --in1/--in2 with --add/--add2)"
        );
    }
    let (fmt, target) =
        io::read_records_with_format(single_input, &args.format, &args.compression)?;
    let add = io::read_records(Some(&args.add), fmt, &args.compression)?;
    let out = spike_single(target, add, args.seed);
    io::write_records(&args.output, fmt, &args.compression, &out)
}

fn spike_single(
    target: Vec<crate::formats::SeqRecord>,
    insert: Vec<crate::formats::SeqRecord>,
    seed: u64,
) -> Vec<crate::formats::SeqRecord> {
    let mut rng = ChaCha20Rng::seed_from_u64(seed);
    let mut inserts = insert
        .into_iter()
        .map(|rec| (rng.random_range(0..=target.len()), rec))
        .collect::<Vec<_>>();
    inserts.sort_unstable_by_key(|(slot, _)| *slot);
    let mut inserts = inserts.into_iter().peekable();
    let mut out = Vec::with_capacity(target.len() + inserts.size_hint().0);
    drain_insert_slot(&mut inserts, 0, &mut out);
    for (i, t) in target.into_iter().enumerate() {
        out.push(t);
        drain_insert_slot(&mut inserts, i + 1, &mut out);
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
    let mut rng = ChaCha20Rng::seed_from_u64(seed);
    let mut inserts = a1
        .into_iter()
        .zip(a2)
        .map(|(x1, x2)| (rng.random_range(0..=t1.len()), x1, x2))
        .collect::<Vec<_>>();
    inserts.sort_unstable_by_key(|(slot, _, _)| *slot);
    let mut inserts = inserts.into_iter().peekable();
    let mut o1 = Vec::with_capacity(t1.len() + inserts.size_hint().0);
    let mut o2 = Vec::with_capacity(t2.len() + inserts.size_hint().0);
    drain_pair_insert_slot(&mut inserts, 0, &mut o1, &mut o2);
    for (i, (x1, x2)) in t1.into_iter().zip(t2.into_iter()).enumerate() {
        o1.push(x1);
        o2.push(x2);
        drain_pair_insert_slot(&mut inserts, i + 1, &mut o1, &mut o2);
    }
    (o1, o2)
}

fn drain_insert_slot(
    inserts: &mut std::iter::Peekable<impl Iterator<Item = (usize, crate::formats::SeqRecord)>>,
    slot: usize,
    out: &mut Vec<crate::formats::SeqRecord>,
) {
    while inserts.peek().is_some_and(|(s, _)| *s == slot) {
        let (_, rec) = inserts.next().expect("peeked insert");
        out.push(rec);
    }
}

fn drain_pair_insert_slot(
    inserts: &mut std::iter::Peekable<
        impl Iterator<Item = (usize, crate::formats::SeqRecord, crate::formats::SeqRecord)>,
    >,
    slot: usize,
    out1: &mut Vec<crate::formats::SeqRecord>,
    out2: &mut Vec<crate::formats::SeqRecord>,
) {
    while inserts.peek().is_some_and(|(s, _, _)| *s == slot) {
        let (_, rec1, rec2) = inserts.next().expect("peeked insert");
        out1.push(rec1);
        out2.push(rec2);
    }
}
