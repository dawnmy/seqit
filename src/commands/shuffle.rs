use anyhow::{Context, Result};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use std::io::Write;

use crate::{cli::ShuffleArgs, formats::SeqFormat, io, pairs, utils};

pub fn run(args: ShuffleArgs) -> Result<()> {
    let mut rng = ChaCha20Rng::seed_from_u64(args.seed);
    let paired_in1 = args.input1.as_deref().or(args.in1.as_deref());
    let paired_in2 = args.input2.as_deref().or(args.in2.as_deref());
    utils::validate_input_mode("shuffle", args.io.input.as_deref(), paired_in1, paired_in2)?;
    if let (Some(in1), Some(in2)) = (paired_in1, paired_in2) {
        let r1 = io::read_records(Some(in1), SeqFormat::Fastq, &args.io.compression)?;
        let r2 = io::read_records(Some(in2), SeqFormat::Fastq, &args.io.compression)?;
        let (r1, r2) = pairs::prepare_paired_records(r1, r2, args.allow_unpaired)?;
        let out2 = args
            .output2
            .as_deref()
            .context("paired shuffle requires --output2")?;
        let mut paired = r1.into_iter().zip(r2).collect::<Vec<_>>();
        use rand::seq::SliceRandom;
        paired.shuffle(&mut rng);
        let w1 = io::open_writer(&args.io.output)?;
        let w1 = io::wrap_compress(w1, &args.io.output, &args.io.compression)?;
        let w2 = io::open_writer(out2)?;
        let w2 = io::wrap_compress(w2, out2, &args.io.compression)?;
        let mut w1 = io::buffered_writer(w1);
        let mut w2 = io::buffered_writer(w2);
        for (a, b) in paired {
            io::write_record(&mut w1, SeqFormat::Fastq, &a)?;
            io::write_record(&mut w2, SeqFormat::Fastq, &b)?;
        }
        w1.flush()?;
        w2.flush()?;
        return Ok(());
    }
    let in_path = args.io.input.as_deref();
    let (fmt, mut recs) =
        io::read_records_with_format(in_path, &args.io.format, &args.io.compression)?;
    use rand::seq::SliceRandom;
    recs.shuffle(&mut rng);
    io::write_records(&args.io.output, fmt, &args.io.compression, &recs)
}
