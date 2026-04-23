use anyhow::{Context, Result};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

use crate::{cli::ShuffleArgs, formats::SeqFormat, io, pairs};

pub fn run(args: ShuffleArgs) -> Result<()> {
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
            .context("paired shuffle requires --output2")?;
        let mut idx: Vec<usize> = (0..r1.len()).collect();
        use rand::seq::SliceRandom;
        idx.shuffle(&mut rng);
        let o1: Vec<_> = idx.iter().map(|i| r1[*i].clone()).collect();
        let o2: Vec<_> = idx.iter().map(|i| r2[*i].clone()).collect();
        io::write_records(&args.io.output, SeqFormat::Fastq, &args.io.compression, &o1)?;
        io::write_records(out2, SeqFormat::Fastq, &args.io.compression, &o2)?;
        return Ok(());
    }
    let in_path = args.io.input.as_deref();
    let (fmt, mut recs) =
        io::read_records_with_format(in_path, &args.io.format, &args.io.compression)?;
    use rand::seq::SliceRandom;
    recs.shuffle(&mut rng);
    io::write_records(&args.io.output, fmt, &args.io.compression, &recs)
}
