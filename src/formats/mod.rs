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
    pub fn len(&self) -> usize {
        self.seq.len()
    }
    pub fn name(&self) -> String {
        match &self.desc {
            Some(d) => format!("{} {}", self.id, d),
            None => self.id.clone(),
        }
    }
}
