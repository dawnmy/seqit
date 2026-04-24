use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};

use anyhow::{bail, Context, Result};
use bio::io::{fasta, fastq};

use crate::cli::{CompressionArg, FormatArg};
use crate::formats::{SeqFormat, SeqRecord};

pub const BUFFER_SIZE: usize = 1024 * 1024;

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

pub fn buffered_reader(reader: Box<dyn Read>) -> BufReader<Box<dyn Read>> {
    BufReader::with_capacity(BUFFER_SIZE, reader)
}

pub fn buffered_writer(writer: Box<dyn Write>) -> BufWriter<Box<dyn Write>> {
    BufWriter::with_capacity(BUFFER_SIZE, writer)
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
    let br = buffered_reader(r);
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
        let r = open_reader(path)?;
        let r = wrap_decompress(r, path, compression)?;
        let mut br = buffered_reader(r);
        let fmt = detect_format_from_bufread(&mut br)?;
        let recs = match fmt {
            SeqFormat::Fasta => read_fasta(br)?,
            SeqFormat::Fastq => read_fastq(br)?,
            _ => bail!("format {:?} not yet implemented for this command", fmt),
        };
        return Ok((fmt, recs));
    }
    let fmt = resolve_seq_format(path, format, compression)?;
    let recs = read_records(path, fmt, compression)?;
    Ok((fmt, recs))
}

pub fn open_seq_reader(
    path: Option<&str>,
    format: &FormatArg,
    compression: &CompressionArg,
) -> Result<(SeqFormat, BufReader<Box<dyn Read>>)> {
    let r = open_reader(path)?;
    let r = wrap_decompress(r, path, compression)?;
    let mut br = buffered_reader(r);

    let fmt = if let Some(fmt) = SeqFormat::from_arg(format) {
        fmt
    } else if let Ok(fmt) = SeqFormat::detect(path) {
        if matches!(fmt, SeqFormat::Fasta | SeqFormat::Fastq) {
            fmt
        } else {
            detect_format_from_bufread(&mut br)?
        }
    } else {
        detect_format_from_bufread(&mut br)?
    };

    Ok((fmt, br))
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
    let mut br = buffered_reader(r);
    detect_format_from_bufread(&mut br)
}

pub(crate) fn detect_format_from_bufread(reader: &mut impl BufRead) -> Result<SeqFormat> {
    let mut first = None;
    loop {
        let buf = reader.fill_buf()?;
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
        reader.consume(consumed);
    }

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
    let mut bw = buffered_writer(w);
    for rec in records {
        write_record(&mut bw, format, rec)?;
    }
    bw.flush()?;
    Ok(())
}

pub fn for_each_record(
    path: Option<&str>,
    format: &FormatArg,
    compression: &CompressionArg,
    mut f: impl FnMut(SeqRecord) -> Result<()>,
) -> Result<SeqFormat> {
    let (fmt, br) = open_seq_reader(path, format, compression)?;
    match fmt {
        SeqFormat::Fasta => {
            let fa = fasta::Reader::new(br);
            for r in fa.records() {
                let r = r?;
                f(SeqRecord {
                    id: r.id().to_string(),
                    desc: r.desc().map(|d| d.to_string()),
                    seq: r.seq().to_vec(),
                    qual: None,
                })?;
            }
        }
        SeqFormat::Fastq => {
            let fq = fastq::Reader::new(br);
            for r in fq.records() {
                let r = r?;
                f(SeqRecord {
                    id: r.id().to_string(),
                    desc: r.desc().map(|d| d.to_string()),
                    seq: r.seq().to_vec(),
                    qual: Some(r.qual().to_vec()),
                })?;
            }
        }
        _ => bail!("format {:?} not yet implemented for this command", fmt),
    }
    Ok(fmt)
}

pub fn write_record(writer: &mut impl Write, format: SeqFormat, rec: &SeqRecord) -> Result<()> {
    match format {
        SeqFormat::Fasta => {
            writer.write_all(b">")?;
            write_header(writer, rec)?;
            writer.write_all(b"\n")?;
            writer.write_all(&rec.seq)?;
            writer.write_all(b"\n")?;
        }
        SeqFormat::Fastq => {
            let Some(q) = &rec.qual else {
                bail!("record '{}' missing qualities for FASTQ output", rec.id);
            };
            writer.write_all(b"@")?;
            write_header(writer, rec)?;
            writer.write_all(b"\n")?;
            writer.write_all(&rec.seq)?;
            writer.write_all(b"\n+\n")?;
            writer.write_all(q)?;
            writer.write_all(b"\n")?;
        }
        _ => bail!("format {:?} not yet implemented for this command", format),
    }
    Ok(())
}

fn write_header(writer: &mut impl Write, rec: &SeqRecord) -> Result<()> {
    writer.write_all(rec.id.as_bytes())?;
    if let Some(desc) = &rec.desc {
        writer.write_all(b" ")?;
        writer.write_all(desc.as_bytes())?;
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
