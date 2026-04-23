use anyhow::{bail, Context, Result};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;

use crate::{cli::SampleArgs, formats::SeqFormat, io, pairs};

pub fn run(args: SampleArgs) -> Result<()> {
    if args.num.is_none() && args.rate.is_none() {
        bail!("one of --num or --rate is required");
    }
    let mut rng = ChaCha20Rng::seed_from_u64(args.seed);

    let paired_in1 = args.input1.as_deref().or(args.in1.as_deref());
    let paired_in2 = args.input2.as_deref().or(args.in2.as_deref());
    if let (Some(in1), Some(in2)) = (paired_in1, paired_in2) {
        let r1 = io::read_records(Some(in1), SeqFormat::Fastq, &args.io.compression)?;
        let r2 = io::read_records(Some(in2), SeqFormat::Fastq, &args.io.compression)?;
        let (r1, r2) = pairs::prepare_paired_records(r1, r2, args.allow_unpaired)?;
        let out2 = args
            .output2
            .as_deref()
            .context("paired sample requires --output2")?;
        if let Some(rate) = args.rate {
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
            return Ok(());
        }
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
        return Ok(());
    }

    let in_path = args.io.input.as_deref();
    let fmt = SeqFormat::from_arg_or_detect(&args.io.format, in_path)?;
    let recs = io::read_records(in_path, fmt, &args.io.compression)?;
    let out = if let Some(rate) = args.rate {
        recs.into_iter()
            .filter(|_| rng.random::<f64>() <= rate)
            .collect()
    } else {
        reservoir(recs, args.num.unwrap(), &mut rng)
    };
    io::write_records(&args.io.output, fmt, &args.io.compression, &out)
}

fn reservoir<T: Clone>(items: Vec<T>, n: usize, rng: &mut ChaCha20Rng) -> Vec<T> {
    if n == 0 {
        return Vec::new();
    }
    let mut out = Vec::new();
    for (i, item) in items.into_iter().enumerate() {
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
