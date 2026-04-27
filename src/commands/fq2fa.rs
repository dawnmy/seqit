use crate::cli::Fq2faArgs;
use crate::formats::SeqFormat;
use crate::io;
use anyhow::{bail, Result};
use bio::io::{fasta, fastq};

pub fn run(args: Fq2faArgs) -> Result<()> {
    let in_path = args.io.input.as_deref();
    let (fmt, reader) = io::open_seq_reader(in_path, &args.io.format, &args.io.compression)?;
    if fmt != SeqFormat::Fastq {
        bail!("fq2fa requires FASTQ input");
    }

    let fq = fastq::Reader::new(reader);
    let w = io::open_writer(&args.io.output)?;
    let w = io::wrap_compress(w, &args.io.output, &args.io.compression)?;
    let mut fa = fasta::Writer::new(io::buffered_writer(w));
    for rec in fq.records() {
        let rec = rec?;
        fa.write(rec.id(), rec.desc(), rec.seq())?;
    }
    Ok(())
}
