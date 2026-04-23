use anyhow::{bail, Result};
use rayon::prelude::*;

use crate::cli::Fq2faArgs;
use crate::formats::SeqFormat;
use crate::io;

pub fn run(args: Fq2faArgs) -> Result<()> {
    let in_path = args.io.input.as_deref();
    let mut recs = if let Some(fmt) = SeqFormat::from_arg(&args.io.format) {
        if fmt != SeqFormat::Fastq {
            bail!("fq2fa requires FASTQ input")
        }
        io::read_records(in_path, SeqFormat::Fastq, &args.io.compression)?
    } else {
        let (fmt, recs) =
            io::read_records_with_format(in_path, &args.io.format, &args.io.compression)?;
        if fmt != SeqFormat::Fastq {
            bail!("fq2fa requires FASTQ input")
        }
        recs
    };
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
