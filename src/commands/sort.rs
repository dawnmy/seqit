use std::{
    cmp::Ordering,
    collections::BinaryHeap,
    fs::{self, File},
    io::{self as stdio, BufReader, BufWriter, ErrorKind, Read, Write},
    path::PathBuf,
    process,
    time::{SystemTime, UNIX_EPOCH},
};

use anyhow::{bail, Context, Result};
use rayon::slice::ParallelSliceMut;

use crate::{
    cli::{SortArgs, SortBy},
    formats::{SeqFormat, SeqRecord},
    io as seqio,
    utils::parse_mem_bytes,
};

pub fn run(args: SortArgs) -> Result<()> {
    let in_path = args.io.input.as_deref();
    let spec = SortSpec {
        by: args.by.clone(),
        numeric: args.numeric,
        reverse: args.reverse,
    };
    let mem_limit = parse_mem_bytes(&args.mem).max(1);
    let tmp_dir = args
        .tmp_dir
        .as_deref()
        .map(PathBuf::from)
        .unwrap_or_else(std::env::temp_dir);
    let mut sorter = ExternalSorter::new(spec, mem_limit, tmp_dir);
    let fmt = seqio::for_each_record(in_path, &args.io.format, &args.io.compression, |rec| {
        sorter.push(rec)
    })?;
    sorter.finish(fmt, &args)
}

#[derive(Clone)]
struct SortSpec {
    by: SortBy,
    numeric: bool,
    reverse: bool,
}

struct ExternalSorter {
    spec: SortSpec,
    mem_limit: usize,
    chunk: Vec<SeqRecord>,
    chunk_bytes: usize,
    chunk_paths: Vec<PathBuf>,
    tmp_dir: PathBuf,
    next_chunk: usize,
    nonce: u128,
}

impl ExternalSorter {
    fn new(spec: SortSpec, mem_limit: usize, tmp_dir: PathBuf) -> Self {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0);
        Self {
            spec,
            mem_limit,
            chunk: Vec::new(),
            chunk_bytes: 0,
            chunk_paths: Vec::new(),
            tmp_dir,
            next_chunk: 0,
            nonce,
        }
    }

    fn push(&mut self, rec: SeqRecord) -> Result<()> {
        self.chunk_bytes = self.chunk_bytes.saturating_add(rec.approx_bytes());
        self.chunk.push(rec);
        if self.chunk_bytes >= self.mem_limit {
            self.flush_chunk()?;
        }
        Ok(())
    }

    fn finish(mut self, fmt: SeqFormat, args: &SortArgs) -> Result<()> {
        if self.chunk_paths.is_empty() {
            sort_records(&mut self.chunk, &self.spec);
            return seqio::write_records(&args.io.output, fmt, &args.io.compression, &self.chunk);
        }
        if !self.chunk.is_empty() {
            self.flush_chunk()?;
        }
        let result = merge_chunks(&self.chunk_paths, fmt, &self.spec, args);
        for path in &self.chunk_paths {
            let _ = fs::remove_file(path);
        }
        result
    }

    fn flush_chunk(&mut self) -> Result<()> {
        if self.chunk.is_empty() {
            return Ok(());
        }
        fs::create_dir_all(&self.tmp_dir)
            .with_context(|| format!("failed to create temp dir '{}'", self.tmp_dir.display()))?;
        sort_records(&mut self.chunk, &self.spec);
        let path = self.next_chunk_path();
        self.next_chunk += 1;
        let file = File::create(&path)
            .with_context(|| format!("failed to create sort chunk '{}'", path.display()))?;
        let mut writer = BufWriter::with_capacity(seqio::BUFFER_SIZE, file);
        for rec in self.chunk.drain(..) {
            write_binary_record(&mut writer, &rec)?;
        }
        writer.flush()?;
        self.chunk_bytes = 0;
        self.chunk_paths.push(path);
        Ok(())
    }

    fn next_chunk_path(&self) -> PathBuf {
        self.tmp_dir.join(format!(
            "seqit-sort-{}-{}-{}.bin",
            process::id(),
            self.nonce,
            self.next_chunk
        ))
    }
}

fn sort_records(records: &mut [SeqRecord], spec: &SortSpec) {
    records.par_sort_unstable_by(|a, b| compare_records(a, b, spec));
}

fn merge_chunks(paths: &[PathBuf], fmt: SeqFormat, spec: &SortSpec, args: &SortArgs) -> Result<()> {
    let mut readers = Vec::with_capacity(paths.len());
    let mut heap = BinaryHeap::new();
    for (chunk_idx, path) in paths.iter().enumerate() {
        let mut reader = ChunkReader::open(path)?;
        if let Some(rec) = reader.next_record()? {
            heap.push(HeapItem::new(rec, chunk_idx, spec));
        }
        readers.push(reader);
    }

    let writer = seqio::open_writer(&args.io.output)?;
    let writer = seqio::wrap_compress(writer, &args.io.output, &args.io.compression)?;
    let mut writer = seqio::buffered_writer(writer);
    while let Some(item) = heap.pop() {
        seqio::write_record(&mut writer, fmt, &item.rec)?;
        if let Some(next) = readers[item.chunk_idx].next_record()? {
            heap.push(HeapItem::new(next, item.chunk_idx, spec));
        }
    }
    writer.flush()?;
    Ok(())
}

struct ChunkReader {
    reader: BufReader<File>,
}

impl ChunkReader {
    fn open(path: &PathBuf) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("failed to open sort chunk '{}'", path.display()))?;
        Ok(Self {
            reader: BufReader::with_capacity(seqio::BUFFER_SIZE, file),
        })
    }

    fn next_record(&mut self) -> Result<Option<SeqRecord>> {
        read_binary_record(&mut self.reader)
    }
}

struct HeapItem {
    key: SortKey,
    rec: SeqRecord,
    chunk_idx: usize,
    reverse: bool,
}

impl HeapItem {
    fn new(rec: SeqRecord, chunk_idx: usize, spec: &SortSpec) -> Self {
        Self {
            key: SortKey::from_record(&rec, spec),
            rec,
            chunk_idx,
            reverse: spec.reverse,
        }
    }
}

impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        let ord = compare_keys(&self.key, &other.key);
        let ord = if self.reverse { ord } else { ord.reverse() };
        ord.then_with(|| other.chunk_idx.cmp(&self.chunk_idx))
    }
}

impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for HeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key && self.chunk_idx == other.chunk_idx
    }
}

impl Eq for HeapItem {}

#[derive(Debug, PartialEq, Eq)]
enum SortKey {
    Len(usize),
    Bytes(Vec<u8>),
    Numeric {
        parsed: Option<u128>,
        bytes: Vec<u8>,
    },
}

impl SortKey {
    fn from_record(rec: &SeqRecord, spec: &SortSpec) -> Self {
        match spec.by {
            SortBy::Len => Self::Len(rec.seq.len()),
            SortBy::Seq => Self::Bytes(rec.seq.clone()),
            SortBy::Id if spec.numeric => {
                let bytes = rec.id.as_bytes().to_vec();
                Self::Numeric {
                    parsed: rec.id.parse::<u128>().ok(),
                    bytes,
                }
            }
            SortBy::Id => Self::Bytes(rec.id.as_bytes().to_vec()),
            SortBy::Name if spec.numeric => {
                let bytes = name_bytes(rec);
                Self::Numeric {
                    parsed: std::str::from_utf8(&bytes)
                        .ok()
                        .and_then(|s| s.parse::<u128>().ok()),
                    bytes,
                }
            }
            SortBy::Name => Self::Bytes(name_bytes(rec)),
        }
    }
}

fn compare_keys(a: &SortKey, b: &SortKey) -> Ordering {
    match (a, b) {
        (SortKey::Len(a), SortKey::Len(b)) => a.cmp(b),
        (SortKey::Bytes(a), SortKey::Bytes(b)) => a.cmp(b),
        (
            SortKey::Numeric {
                parsed: Some(a), ..
            },
            SortKey::Numeric {
                parsed: Some(b), ..
            },
        ) => a.cmp(b),
        (SortKey::Numeric { bytes: a, .. }, SortKey::Numeric { bytes: b, .. }) => a.cmp(b),
        _ => Ordering::Equal,
    }
}

fn compare_records(a: &SeqRecord, b: &SeqRecord, spec: &SortSpec) -> Ordering {
    let ord = match spec.by {
        SortBy::Len => a.seq.len().cmp(&b.seq.len()),
        SortBy::Seq => a.seq.cmp(&b.seq),
        SortBy::Id if spec.numeric => compare_numeric_str(&a.id, &b.id),
        SortBy::Id => a.id.cmp(&b.id),
        SortBy::Name if spec.numeric => compare_numeric_bytes(name_bytes(a), name_bytes(b)),
        SortBy::Name => name_iter(a).cmp(name_iter(b)),
    };
    if spec.reverse {
        ord.reverse()
    } else {
        ord
    }
}

fn compare_numeric_str(a: &str, b: &str) -> Ordering {
    match (a.parse::<u128>(), b.parse::<u128>()) {
        (Ok(a), Ok(b)) => a.cmp(&b),
        _ => a.as_bytes().cmp(b.as_bytes()),
    }
}

fn compare_numeric_bytes(a: Vec<u8>, b: Vec<u8>) -> Ordering {
    match (
        std::str::from_utf8(&a)
            .ok()
            .and_then(|s| s.parse::<u128>().ok()),
        std::str::from_utf8(&b)
            .ok()
            .and_then(|s| s.parse::<u128>().ok()),
    ) {
        (Some(a), Some(b)) => a.cmp(&b),
        _ => a.cmp(&b),
    }
}

fn name_iter(rec: &SeqRecord) -> impl Iterator<Item = u8> + '_ {
    rec.id.bytes().chain(
        rec.desc
            .iter()
            .flat_map(|desc| std::iter::once(b' ').chain(desc.bytes())),
    )
}

fn name_bytes(rec: &SeqRecord) -> Vec<u8> {
    let mut out = Vec::with_capacity(rec.id.len() + rec.desc.as_ref().map_or(0, |d| d.len() + 1));
    out.extend_from_slice(rec.id.as_bytes());
    if let Some(desc) = &rec.desc {
        out.push(b' ');
        out.extend_from_slice(desc.as_bytes());
    }
    out
}

const NONE_LEN: u64 = u64::MAX;

fn write_binary_record(writer: &mut impl Write, rec: &SeqRecord) -> Result<()> {
    write_bytes(writer, rec.id.as_bytes())?;
    match &rec.desc {
        Some(desc) => write_bytes(writer, desc.as_bytes())?,
        None => write_len(writer, NONE_LEN)?,
    }
    write_bytes(writer, &rec.seq)?;
    match &rec.qual {
        Some(qual) => write_bytes(writer, qual)?,
        None => write_len(writer, NONE_LEN)?,
    }
    Ok(())
}

fn read_binary_record(reader: &mut impl Read) -> Result<Option<SeqRecord>> {
    let Some(id_len) = read_len_or_eof(reader)? else {
        return Ok(None);
    };
    let id = String::from_utf8(read_exact_len(reader, id_len)?)
        .context("sort chunk contains invalid UTF-8 record id")?;
    let desc_len = read_len(reader)?;
    let desc = if desc_len == NONE_LEN {
        None
    } else {
        Some(
            String::from_utf8(read_exact_len(reader, desc_len)?)
                .context("sort chunk contains invalid UTF-8 record description")?,
        )
    };
    let seq_len = read_len(reader)?;
    let seq = read_exact_len(reader, seq_len)?;
    let qual_len = read_len(reader)?;
    let qual = if qual_len == NONE_LEN {
        None
    } else {
        Some(read_exact_len(reader, qual_len)?)
    };
    Ok(Some(SeqRecord {
        id,
        desc,
        seq,
        qual,
    }))
}

fn write_bytes(writer: &mut impl Write, bytes: &[u8]) -> Result<()> {
    write_len(writer, bytes.len() as u64)?;
    writer.write_all(bytes)?;
    Ok(())
}

fn write_len(writer: &mut impl Write, len: u64) -> Result<()> {
    writer.write_all(&len.to_le_bytes())?;
    Ok(())
}

fn read_len_or_eof(reader: &mut impl Read) -> Result<Option<u64>> {
    let mut buf = [0u8; 8];
    match reader.read_exact(&mut buf) {
        Ok(()) => Ok(Some(u64::from_le_bytes(buf))),
        Err(err) if err.kind() == ErrorKind::UnexpectedEof => Ok(None),
        Err(err) => Err(err.into()),
    }
}

fn read_len(reader: &mut impl Read) -> Result<u64> {
    read_len_or_eof(reader)?.ok_or_else(|| {
        stdio::Error::new(
            ErrorKind::UnexpectedEof,
            "truncated seqit sort chunk record",
        )
        .into()
    })
}

fn read_exact_len(reader: &mut impl Read, len: u64) -> Result<Vec<u8>> {
    if len > usize::MAX as u64 {
        bail!("sort chunk record is too large for this platform");
    }
    let mut out = vec![0u8; len as usize];
    reader.read_exact(&mut out)?;
    Ok(out)
}
