use anyhow::{bail, Result};
use bio::io::{fasta, fastq};
use std::io::BufReader;

use crate::cli::{FormatArg, Fq2faArgs};
use crate::formats::SeqFormat;
use crate::io;

pub fn run(args: Fq2faArgs) -> Result<()> {
    let in_path = args.io.input.as_deref();
    let fmt = if matches!(args.io.format, FormatArg::Auto) && matches!(in_path, None | Some("-")) {
        let (fmt, recs) =
            io::read_records_with_format(in_path, &args.io.format, &args.io.compression)?;
        if fmt != SeqFormat::Fastq {
            bail!("fq2fa requires FASTQ input");
        }
        let mut out = recs;
        for rec in &mut out {
            rec.qual = None;
        }
        return io::write_records(
            &args.io.output,
            SeqFormat::Fasta,
            &args.io.compression,
            &out,
        );
    } else {
        io::resolve_seq_format(in_path, &args.io.format, &args.io.compression)?
    };
    if fmt != SeqFormat::Fastq {
        bail!("fq2fa requires FASTQ input");
    }

    let r = io::open_reader(in_path)?;
    let r = io::wrap_decompress(r, in_path, &args.io.compression)?;
    let fq = fastq::Reader::new(BufReader::new(r));
    let w = io::open_writer(&args.io.output)?;
    let w = io::wrap_compress(w, &args.io.output, &args.io.compression)?;
    let mut fa = fasta::Writer::new(w);
    for rec in fq.records() {
        let rec = rec?;
        fa.write(rec.id(), rec.desc(), rec.seq())?;
    }
    Ok(())
}
