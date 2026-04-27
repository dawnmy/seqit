use crate::cli::FormatArg;
use anyhow::{bail, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqFormat {
    Fasta,
    Fastq,
    Sam,
    Bam,
    Cram,
}

impl SeqFormat {
    pub fn from_arg(arg: &FormatArg) -> Option<Self> {
        match arg {
            FormatArg::Auto => None,
            FormatArg::Fasta => Some(Self::Fasta),
            FormatArg::Fastq => Some(Self::Fastq),
            FormatArg::Sam => Some(Self::Sam),
            FormatArg::Bam => Some(Self::Bam),
            FormatArg::Cram => Some(Self::Cram),
        }
    }

    pub fn from_arg_or_detect(arg: &FormatArg, path: Option<&str>) -> Result<Self> {
        if let Some(fmt) = Self::from_arg(arg) {
            Ok(fmt)
        } else {
            Self::detect(path)
        }
    }

    pub fn detect(path: Option<&str>) -> Result<Self> {
        let p = path.unwrap_or("-");
        let stripped = p
            .strip_suffix(".gz")
            .or_else(|| p.strip_suffix(".xz"))
            .unwrap_or(p);
        if stripped.ends_with(".fa") || stripped.ends_with(".fasta") || stripped.ends_with(".fna") {
            Ok(Self::Fasta)
        } else if stripped.ends_with(".fq") || stripped.ends_with(".fastq") {
            Ok(Self::Fastq)
        } else if stripped.ends_with(".sam") {
            Ok(Self::Sam)
        } else if stripped.ends_with(".bam") {
            Ok(Self::Bam)
        } else if stripped.ends_with(".cram") {
            Ok(Self::Cram)
        } else {
            bail!("Unable to auto-detect format for '{p}'. Use --format.")
        }
    }
}

#[derive(Debug, Clone)]
pub struct SeqRecord {
    pub id: String,
    pub desc: Option<String>,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>,
}

impl SeqRecord {
    pub fn approx_bytes(&self) -> usize {
        const RECORD_OVERHEAD: usize = 96;
        RECORD_OVERHEAD
            + self.id.len()
            + self.desc.as_ref().map_or(0, |d| d.len())
            + self.seq.len()
            + self.qual.as_ref().map_or(0, |q| q.len())
    }
}
