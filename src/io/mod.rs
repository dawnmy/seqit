use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};

use anyhow::{bail, Context, Result};
use bio::io::{fasta, fastq};

use crate::cli::CompressionArg;
use crate::formats::{SeqFormat, SeqRecord};

pub fn open_reader(path: Option<&str>) -> Result<Box<dyn Read>> {
    let reader: Box<dyn Read> = match path {
        None | Some("-") => Box::new(io::stdin()),
        Some(p) => Box::new(File::open(p).with_context(|| format!("failed to open input '{p}'"))?),
    };
    Ok(reader)
}

pub fn open_writer(path: &str) -> Result<Box<dyn Write>> {
    let writer: Box<dyn Write> = if path == "-" {
        Box::new(io::stdout())
    } else {
        Box::new(File::create(path).with_context(|| format!("failed to create output '{path}'"))?)
    };
    Ok(writer)
}

pub fn wrap_decompress(
    reader: Box<dyn Read>,
    path: Option<&str>,
    mode: &CompressionArg,
) -> Result<Box<dyn Read>> {
    match detect_compression(path, mode) {
        CompressionArg::Gz => Ok(Box::new(flate2::read::MultiGzDecoder::new(reader))),
        CompressionArg::Xz => Ok(Box::new(xz2::read::XzDecoder::new(reader))),
        _ => Ok(reader),
    }
}

pub fn wrap_compress(
    writer: Box<dyn Write>,
    path: &str,
    mode: &CompressionArg,
) -> Result<Box<dyn Write>> {
    match detect_compression(Some(path), mode) {
        CompressionArg::Gz => Ok(Box::new(flate2::write::GzEncoder::new(
            writer,
            flate2::Compression::default(),
        ))),
        CompressionArg::Xz => Ok(Box::new(xz2::write::XzEncoder::new(writer, 6))),
        _ => Ok(writer),
    }
}

fn detect_compression(path: Option<&str>, mode: &CompressionArg) -> CompressionArg {
    match mode {
        CompressionArg::Auto => {
            let p = path.unwrap_or("-");
            if p.ends_with(".gz") {
                CompressionArg::Gz
            } else if p.ends_with(".xz") {
                CompressionArg::Xz
            } else {
                CompressionArg::None
            }
        }
        other => other.clone(),
    }
}

pub fn read_records(
    path: Option<&str>,
    format: SeqFormat,
    compression: &CompressionArg,
) -> Result<Vec<SeqRecord>> {
    let r = open_reader(path)?;
    let r = wrap_decompress(r, path, compression)?;
    let br = BufReader::new(r);
    match format {
        SeqFormat::Fasta => read_fasta(br),
        SeqFormat::Fastq => read_fastq(br),
        _ => bail!("format {:?} not yet implemented for this command", format),
    }
}

pub fn write_records(
    path: &str,
    format: SeqFormat,
    compression: &CompressionArg,
    records: &[SeqRecord],
) -> Result<()> {
    let w = open_writer(path)?;
    let w = wrap_compress(w, path, compression)?;
    let bw = BufWriter::new(w);
    match format {
        SeqFormat::Fasta => write_fasta(bw, records),
        SeqFormat::Fastq => write_fastq(bw, records),
        _ => bail!("format {:?} not yet implemented for this command", format),
    }
}

pub fn stream_records<F>(
    path: Option<&str>,
    format: SeqFormat,
    compression: &CompressionArg,
    mut f: F,
) -> Result<()>
where
    F: FnMut(SeqRecord) -> Result<()>,
{
    for rec in read_records(path, format, compression)? {
        f(rec)?;
    }
    Ok(())
}

fn read_fasta(reader: impl BufRead) -> Result<Vec<SeqRecord>> {
    let mut out = Vec::new();
    let fa = fasta::Reader::new(reader);
    for r in fa.records() {
        let r = r?;
        out.push(SeqRecord {
            id: r.id().to_string(),
            desc: r.desc().map(|d| d.to_string()),
            seq: r.seq().to_vec(),
            qual: None,
        });
    }
    Ok(out)
}

fn read_fastq(reader: impl BufRead) -> Result<Vec<SeqRecord>> {
    let mut out = Vec::new();
    let fq = fastq::Reader::new(reader);
    for r in fq.records() {
        let r = r?;
        out.push(SeqRecord {
            id: r.id().to_string(),
            desc: r.desc().map(|d| d.to_string()),
            seq: r.seq().to_vec(),
            qual: Some(r.qual().to_vec()),
        });
    }
    Ok(out)
}

fn write_fasta(writer: impl Write, records: &[SeqRecord]) -> Result<()> {
    let mut w = fasta::Writer::new(writer);
    for r in records {
        w.write(&r.id, r.desc.as_deref(), &r.seq)?;
    }
    Ok(())
}

fn write_fastq(writer: impl Write, records: &[SeqRecord]) -> Result<()> {
    let mut w = fastq::Writer::new(writer);
    for r in records {
        let Some(q) = &r.qual else {
            bail!("record '{}' missing qualities for FASTQ output", r.id);
        };
        w.write(&r.id, r.desc.as_deref(), &r.seq, q)?;
    }
    Ok(())
}
