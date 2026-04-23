use anyhow::{bail, Result};
use rayon::prelude::*;

use crate::cli::Fq2faArgs;
use crate::formats::SeqFormat;
use crate::io;

pub fn run(args: Fq2faArgs) -> Result<()> {
    let in_path = args.io.input.as_deref();
    let fmt = SeqFormat::from_arg_or_detect(&args.io.format, in_path)?;
    if fmt != SeqFormat::Fastq {
        bail!("fq2fa requires FASTQ input")
    }
    let mut recs = io::read_records(in_path, SeqFormat::Fastq, &args.io.compression)?;
    recs.par_iter_mut().for_each(|r| {
        r.qual = None;
    });
    io::write_records(
        &args.io.output,
        SeqFormat::Fasta,
        &args.io.compression,
        &recs,
    )
}
