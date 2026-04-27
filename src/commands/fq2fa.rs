use crate::cli::Fq2faArgs;
use crate::formats::SeqFormat;
use crate::io;
use anyhow::{bail, Result};
use std::io::Write;

pub fn run(args: Fq2faArgs) -> Result<()> {
    let in_path = args.io.input.as_deref();
    let (fmt, reader) = io::open_seq_reader(in_path, &args.io.format, &args.io.compression)?;
    if fmt != SeqFormat::Fastq {
        bail!("fq2fa requires FASTQ input");
    }

    let w = io::open_writer(&args.io.output)?;
    let w = io::wrap_compress(w, &args.io.output, &args.io.compression)?;
    let mut out = io::buffered_writer(w);
    io::for_each_record_from_reader(reader, fmt, |rec| {
        out.write_all(b">")?;
        io::write_header_bytes(&mut out, rec.id, rec.desc)?;
        out.write_all(b"\n")?;
        out.write_all(&rec.seq)?;
        out.write_all(b"\n")?;
        Ok(())
    })?;
    out.flush()?;
    Ok(())
}
