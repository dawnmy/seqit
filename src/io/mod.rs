use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Cursor, Read, Write};

use anyhow::{bail, Context, Result};
use bio::io::{fasta, fastq};

use crate::cli::{CompressionArg, FormatArg};
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

pub fn read_records_with_format(
    path: Option<&str>,
    format: &FormatArg,
    compression: &CompressionArg,
) -> Result<(SeqFormat, Vec<SeqRecord>)> {
    if matches!(format, FormatArg::Auto) && matches!(path, None | Some("-")) {
        return read_records_detect_from_content(path, compression);
    }
    let fmt = resolve_seq_format(path, format, compression)?;
    let recs = read_records(path, fmt, compression)?;
    Ok((fmt, recs))
}

pub fn resolve_seq_format(
    path: Option<&str>,
    format: &FormatArg,
    compression: &CompressionArg,
) -> Result<SeqFormat> {
    if let Some(fmt) = SeqFormat::from_arg(format) {
        return Ok(fmt);
    }

    if let Ok(fmt) = SeqFormat::detect(path) {
        if matches!(fmt, SeqFormat::Fasta | SeqFormat::Fastq) {
            return Ok(fmt);
        }
    }

    detect_format_from_stream(path, compression)
}

fn detect_format_from_stream(
    path: Option<&str>,
    compression: &CompressionArg,
) -> Result<SeqFormat> {
    let r = open_reader(path)?;
    let r = wrap_decompress(r, path, compression)?;
    let mut br = BufReader::new(r);
    let mut first = None;
    loop {
        let buf = br.fill_buf()?;
        if buf.is_empty() {
            break;
        }
        if let Some(b) = buf
            .iter()
            .copied()
            .find(|c| !matches!(c, b' ' | b'\n' | b'\r' | b'\t'))
        {
            first = Some(b);
            break;
        }
        let consumed = buf.len();
        br.consume(consumed);
    }
    match first {
        Some(b'>') => Ok(SeqFormat::Fasta),
        Some(b'@') => Ok(SeqFormat::Fastq),
        _ => bail!("Unable to auto-detect format from input content. Use --format."),
    }
}

fn read_records_detect_from_content(
    path: Option<&str>,
    compression: &CompressionArg,
) -> Result<(SeqFormat, Vec<SeqRecord>)> {
    let r = open_reader(path)?;
    let mut r = wrap_decompress(r, path, compression)?;
    let mut buf = Vec::new();
    r.read_to_end(&mut buf)?;

    let fmt = detect_format_from_content(&buf)?;
    let br = BufReader::new(Cursor::new(buf));
    let recs = match fmt {
        SeqFormat::Fasta => read_fasta(br)?,
        SeqFormat::Fastq => read_fastq(br)?,
        _ => bail!("format {:?} not yet implemented for this command", fmt),
    };
    Ok((fmt, recs))
}

fn detect_format_from_content(buf: &[u8]) -> Result<SeqFormat> {
    let first = buf
        .iter()
        .copied()
        .find(|b| !matches!(b, b' ' | b'\n' | b'\r' | b'\t'));
    match first {
        Some(b'>') => Ok(SeqFormat::Fasta),
        Some(b'@') => Ok(SeqFormat::Fastq),
        _ => bail!("Unable to auto-detect format from input content. Use --format."),
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
