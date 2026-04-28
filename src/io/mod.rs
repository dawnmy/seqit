use std::borrow::Cow;
use std::fs::File;
use std::io::{self, BufRead, BufWriter, Read, Write};
use std::str;

use anyhow::{bail, Context, Result};

use crate::cli::{CompressionArg, FormatArg};
use crate::formats::{SeqFormat, SeqRecord};

mod aligned;
mod fastx;

use aligned::AlignedBufReader;

pub const BUFFER_SIZE: usize = 1024 * 1024;
const BUFFER_ALIGNMENT: usize = 64;
const PARALLEL_GZIP_LEVEL: u32 = 3;

pub type DynReader = Box<dyn Read + Send>;
pub type DynBufReader = AlignedBufReader<DynReader>;
pub type DynWriter = Box<dyn Write>;
pub type DynSendWriter = Box<dyn Write + Send>;

pub fn open_reader(path: Option<&str>) -> Result<DynReader> {
    let reader: DynReader = match path {
        None | Some("-") => Box::new(io::stdin()),
        Some(p) => Box::new(File::open(p).with_context(|| format!("failed to open input '{p}'"))?),
    };
    Ok(reader)
}

pub fn open_writer(path: &str) -> Result<DynSendWriter> {
    let writer: DynSendWriter = if path == "-" {
        Box::new(io::stdout())
    } else {
        Box::new(File::create(path).with_context(|| format!("failed to create output '{path}'"))?)
    };
    Ok(writer)
}

pub fn buffered_reader(reader: DynReader) -> DynBufReader {
    AlignedBufReader::with_capacity_alignment(BUFFER_SIZE, BUFFER_ALIGNMENT, reader)
}

pub fn buffered_writer(writer: DynWriter) -> BufWriter<DynWriter> {
    BufWriter::with_capacity(BUFFER_SIZE, writer)
}

pub fn wrap_decompress(
    reader: DynReader,
    path: Option<&str>,
    mode: &CompressionArg,
) -> Result<DynReader> {
    match detect_compression(path, mode) {
        CompressionArg::Gz => Ok(Box::new(flate2::bufread::MultiGzDecoder::new(
            buffered_reader(reader),
        ))),
        CompressionArg::Xz => Ok(Box::new(xz2::bufread::XzDecoder::new_multi_decoder(
            buffered_reader(reader),
        ))),
        _ => Ok(reader),
    }
}

pub fn wrap_compress(
    writer: DynSendWriter,
    path: &str,
    mode: &CompressionArg,
) -> Result<DynWriter> {
    wrap_compress_for_streams(writer, path, mode, 1)
}

pub fn wrap_compress_for_streams(
    writer: DynSendWriter,
    path: &str,
    mode: &CompressionArg,
    concurrent_streams: usize,
) -> Result<DynWriter> {
    match detect_compression(Some(path), mode) {
        CompressionArg::Gz => Ok(wrap_gzip_writer(writer, concurrent_streams)),
        CompressionArg::Xz => wrap_xz_writer(writer, concurrent_streams),
        _ => Ok(writer),
    }
}

fn wrap_gzip_writer(writer: DynSendWriter, concurrent_streams: usize) -> DynWriter {
    let threads = compression_threads_for_streams(concurrent_streams);
    if threads <= 1 {
        return Box::new(flate2::write::GzEncoder::new(
            writer,
            flate2::Compression::default(),
        ));
    }
    Box::new(ParallelGzipWriter::new(writer, threads))
}

fn compression_threads_for_streams(concurrent_streams: usize) -> usize {
    let streams = concurrent_streams.max(1);
    (rayon::current_num_threads() / streams).max(1)
}

fn wrap_xz_writer(writer: DynSendWriter, concurrent_streams: usize) -> Result<DynWriter> {
    let threads = compression_threads_for_streams(concurrent_streams);
    if threads <= 1 {
        return Ok(Box::new(xz2::write::XzEncoder::new(writer, 6)));
    }

    let mut builder = xz2::stream::MtStreamBuilder::new();
    builder
        .threads(threads.min(u32::MAX as usize) as u32)
        .preset(6);
    let stream = builder.encoder()?;
    Ok(Box::new(xz2::write::XzEncoder::new_stream(writer, stream)))
}

struct ParallelGzipWriter {
    inner: Option<Box<dyn gzp::ZWriter<DynSendWriter>>>,
}

impl ParallelGzipWriter {
    fn new(writer: DynSendWriter, threads: usize) -> Self {
        let inner = gzp::ZBuilder::<gzp::deflate::Gzip, _>::new()
            .num_threads(threads)
            .compression_level(flate2::Compression::new(PARALLEL_GZIP_LEVEL))
            .buffer_size(gzp::BUFSIZE)
            .from_writer(writer);
        Self { inner: Some(inner) }
    }
}

impl Write for ParallelGzipWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.inner
            .as_mut()
            .expect("gzip writer used after finish")
            .write(buf)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.inner
            .as_mut()
            .expect("gzip writer used after finish")
            .flush()
    }
}

impl Drop for ParallelGzipWriter {
    fn drop(&mut self) {
        if let Some(mut inner) = self.inner.take() {
            let _ = inner.finish();
        }
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
    read_records_from_reader(br, format)
}

pub fn read_record_pair_parallel(
    in1: &str,
    in2: &str,
    format: SeqFormat,
    compression: &CompressionArg,
) -> Result<(Vec<SeqRecord>, Vec<SeqRecord>)> {
    if can_parallel_read(in1, in2) {
        let (r1, r2) = rayon::join(
            || read_records(Some(in1), format, compression),
            || read_records(Some(in2), format, compression),
        );
        Ok((r1?, r2?))
    } else {
        let r1 = read_records(Some(in1), format, compression)?;
        let r2 = read_records(Some(in2), format, compression)?;
        Ok((r1, r2))
    }
}

fn can_parallel_read(in1: &str, in2: &str) -> bool {
    in1 != "-" && in2 != "-"
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
        let recs = read_records_from_reader(br, fmt)?;
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
) -> Result<(SeqFormat, DynBufReader)> {
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
    write_records_for_streams(path, format, compression, records, 1)
}

fn write_records_for_streams(
    path: &str,
    format: SeqFormat,
    compression: &CompressionArg,
    records: &[SeqRecord],
    concurrent_streams: usize,
) -> Result<()> {
    let w = open_writer(path)?;
    let w = wrap_compress_for_streams(w, path, compression, concurrent_streams)?;
    let mut bw = buffered_writer(w);
    for rec in records {
        write_record(&mut bw, format, rec)?;
    }
    bw.flush()?;
    Ok(())
}

pub fn write_record_pair_parallel(
    out1: &str,
    out2: &str,
    format: SeqFormat,
    compression: &CompressionArg,
    records1: &[SeqRecord],
    records2: &[SeqRecord],
) -> Result<()> {
    if can_parallel_write(out1, out2) {
        let (left, right) = rayon::join(
            || write_records_for_streams(out1, format, compression, records1, 2),
            || write_records_for_streams(out2, format, compression, records2, 2),
        );
        left?;
        right?;
        Ok(())
    } else {
        write_records(out1, format, compression, records1)?;
        write_records(out2, format, compression, records2)
    }
}

fn can_parallel_write(out1: &str, out2: &str) -> bool {
    out1 != "-" && out2 != "-" && out1 != out2
}

pub fn for_each_record(
    path: Option<&str>,
    format: &FormatArg,
    compression: &CompressionArg,
    mut f: impl FnMut(SeqRecord) -> Result<()>,
) -> Result<SeqFormat> {
    for_each_record_ref(path, format, compression, |rec| f(rec.to_owned_record()?))
}

pub fn for_each_record_ref<F>(
    path: Option<&str>,
    format: &FormatArg,
    compression: &CompressionArg,
    f: F,
) -> Result<SeqFormat>
where
    F: for<'a> FnMut(SeqRecordRef<'a>) -> Result<()>,
{
    let (fmt, br) = open_seq_reader(path, format, compression)?;
    for_each_record_from_reader(br, fmt, f)?;
    Ok(fmt)
}

pub(crate) fn for_each_record_from_reader<F>(
    br: DynBufReader,
    fmt: SeqFormat,
    mut f: F,
) -> Result<()>
where
    F: for<'a> FnMut(SeqRecordRef<'a>) -> Result<()>,
{
    for_each_record_from_reader_until(br, fmt, |rec| {
        f(rec)?;
        Ok(true)
    })
}

pub(crate) fn for_each_record_from_reader_until<F>(
    br: DynBufReader,
    fmt: SeqFormat,
    mut f: F,
) -> Result<()>
where
    F: for<'a> FnMut(SeqRecordRef<'a>) -> Result<bool>,
{
    match fmt {
        SeqFormat::Fasta | SeqFormat::Fastq => {
            let mut reader = fastx::Reader::from_reader(br);
            while let Some(record) = reader.next()? {
                if fmt == SeqFormat::Fastq && record.qual.is_none() {
                    bail!("expected FASTQ input but found FASTA records");
                }
                if fmt == SeqFormat::Fasta && record.qual.is_some() {
                    bail!("expected FASTA input but found FASTQ records");
                }
                let rec = SeqRecordRef {
                    id: record.id,
                    desc: record.desc,
                    seq: Cow::Borrowed(record.seq),
                    qual: record.qual,
                };
                if !f(rec)? {
                    break;
                }
            }
        }
        _ => bail!("format {:?} not yet implemented for this command", fmt),
    }
    Ok(())
}

pub(crate) enum FastqPair<'a, 'b> {
    Both(SeqRecordRef<'a>, SeqRecordRef<'b>),
    Left(SeqRecordRef<'a>),
    Right(SeqRecordRef<'b>),
}

pub(crate) fn for_each_fastq_pair<F>(
    in1: &str,
    in2: &str,
    compression: &CompressionArg,
    mut f: F,
) -> Result<()>
where
    F: for<'a, 'b> FnMut(FastqPair<'a, 'b>) -> Result<bool>,
{
    let r1 = open_reader(Some(in1))?;
    let r1 = wrap_decompress(r1, Some(in1), compression)?;
    let r2 = open_reader(Some(in2))?;
    let r2 = wrap_decompress(r2, Some(in2), compression)?;
    let mut r1 = fastx::Reader::from_reader(buffered_reader(r1));
    let mut r2 = fastx::Reader::from_reader(buffered_reader(r2));

    loop {
        let a = r1.next()?;
        let b = r2.next()?;
        match (a, b) {
            (Some(a), Some(b)) => {
                let a = fastq_record_ref(a)?;
                let b = fastq_record_ref(b)?;
                if !f(FastqPair::Both(a, b))? {
                    break;
                }
            }
            (Some(a), None) => {
                let a = fastq_record_ref(a)?;
                if !f(FastqPair::Left(a))? {
                    break;
                }
            }
            (None, Some(b)) => {
                let b = fastq_record_ref(b)?;
                if !f(FastqPair::Right(b))? {
                    break;
                }
            }
            (None, None) => break,
        }
    }
    Ok(())
}

fn fastq_record_ref(record: fastx::Record<'_>) -> Result<SeqRecordRef<'_>> {
    let Some(qual) = record.qual else {
        bail!("expected FASTQ input but found FASTA records");
    };
    Ok(SeqRecordRef {
        id: record.id,
        desc: record.desc,
        seq: Cow::Borrowed(record.seq),
        qual: Some(qual),
    })
}

pub struct SeqRecordRef<'a> {
    pub id: &'a [u8],
    pub desc: Option<&'a [u8]>,
    pub seq: Cow<'a, [u8]>,
    pub qual: Option<&'a [u8]>,
}

impl SeqRecordRef<'_> {
    pub fn id_str(&self) -> Result<&str> {
        bytes_to_str(self.id, "record id")
    }

    pub fn desc_str(&self) -> Result<Option<&str>> {
        self.desc
            .map(|desc| bytes_to_str(desc, "record description"))
            .transpose()
    }

    pub fn to_owned_record(&self) -> Result<SeqRecord> {
        Ok(SeqRecord {
            id: bytes_to_str(self.id, "record id")?.to_string(),
            desc: self
                .desc
                .map(|desc| bytes_to_str(desc, "record description").map(str::to_string))
                .transpose()?,
            seq: self.seq.to_vec(),
            qual: self.qual.map(<[u8]>::to_vec),
        })
    }
}

pub fn bytes_to_str<'a>(bytes: &'a [u8], field: &str) -> Result<&'a str> {
    str::from_utf8(bytes).with_context(|| format!("{field} is not valid UTF-8"))
}

fn read_records_from_reader(reader: DynBufReader, format: SeqFormat) -> Result<Vec<SeqRecord>> {
    let mut out = Vec::new();
    for_each_record_from_reader(reader, format, |rec| {
        out.push(rec.to_owned_record()?);
        Ok(())
    })?;
    Ok(out)
}

pub fn write_record(writer: &mut impl Write, format: SeqFormat, rec: &SeqRecord) -> Result<()> {
    match format {
        SeqFormat::Fasta => {
            writer.write_all(b">")?;
            write_owned_header(writer, rec)?;
            writer.write_all(b"\n")?;
            writer.write_all(&rec.seq)?;
            writer.write_all(b"\n")?;
        }
        SeqFormat::Fastq => {
            let Some(q) = &rec.qual else {
                bail!("record '{}' missing qualities for FASTQ output", rec.id);
            };
            writer.write_all(b"@")?;
            write_owned_header(writer, rec)?;
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

pub fn write_record_ref(
    writer: &mut impl Write,
    format: SeqFormat,
    rec: &SeqRecordRef<'_>,
) -> Result<()> {
    match format {
        SeqFormat::Fasta => {
            writer.write_all(b">")?;
            write_header_bytes(writer, rec.id, rec.desc)?;
            writer.write_all(b"\n")?;
            writer.write_all(&rec.seq)?;
            writer.write_all(b"\n")?;
        }
        SeqFormat::Fastq => {
            let Some(q) = rec.qual else {
                bail!(
                    "record '{}' missing qualities for FASTQ output",
                    rec.id_str()?
                );
            };
            writer.write_all(b"@")?;
            write_header_bytes(writer, rec.id, rec.desc)?;
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

fn write_owned_header(writer: &mut impl Write, rec: &SeqRecord) -> Result<()> {
    writer.write_all(rec.id.as_bytes())?;
    if let Some(desc) = &rec.desc {
        writer.write_all(b" ")?;
        writer.write_all(desc.as_bytes())?;
    }
    Ok(())
}

pub fn write_header_bytes(writer: &mut impl Write, id: &[u8], desc: Option<&[u8]>) -> Result<()> {
    writer.write_all(id)?;
    if let Some(desc) = desc {
        writer.write_all(b" ")?;
        writer.write_all(desc)?;
    }
    Ok(())
}
