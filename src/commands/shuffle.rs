use anyhow::{Context, Result};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

use crate::{cli::ShuffleArgs, formats::SeqFormat, io, pairs, utils};

pub fn run(args: ShuffleArgs) -> Result<()> {
    let mut rng = ChaCha20Rng::seed_from_u64(args.seed);
    let paired_in1 = args.input1.as_deref().or(args.in1.as_deref());
    let paired_in2 = args.input2.as_deref().or(args.in2.as_deref());
    utils::validate_input_mode("shuffle", args.io.input.as_deref(), paired_in1, paired_in2)?;
    if let (Some(in1), Some(in2)) = (paired_in1, paired_in2) {
        let (r1, r2) =
            io::read_record_pair_parallel(in1, in2, SeqFormat::Fastq, &args.io.compression)?;
        let (r1, r2) = pairs::prepare_paired_records(r1, r2, args.allow_unpaired)?;
        let out2 = args
            .output2
            .as_deref()
            .context("paired shuffle requires --output2")?;
        let mut paired = r1.into_iter().zip(r2).collect::<Vec<_>>();
        use rand::seq::SliceRandom;
        paired.shuffle(&mut rng);
        let (o1, o2): (Vec<_>, Vec<_>) = paired.into_iter().unzip();
        io::write_record_pair_parallel(
            &args.io.output,
            out2,
            SeqFormat::Fastq,
            &args.io.compression,
            &o1,
            &o2,
        )?;
        return Ok(());
    }
    let in_path = args.io.input.as_deref();
    let (fmt, mut recs) =
        io::read_records_with_format(in_path, &args.io.format, &args.io.compression)?;
    use rand::seq::SliceRandom;
    recs.shuffle(&mut rng);
    io::write_records(&args.io.output, fmt, &args.io.compression, &recs)
}
