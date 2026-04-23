use ahash::AHasher;
use anyhow::{bail, Result};
use std::hash::{Hash, Hasher};

pub fn hash64<T: Hash>(value: &T) -> u64 {
    let mut h = AHasher::default();
    value.hash(&mut h);
    h.finish()
}

pub fn parse_mem_bytes(mem: &str) -> usize {
    let m = mem.trim().to_ascii_uppercase();
    if let Some(v) = m.strip_suffix('G') {
        return v.parse::<usize>().unwrap_or(1) * 1024 * 1024 * 1024;
    }
    if let Some(v) = m.strip_suffix('M') {
        return v.parse::<usize>().unwrap_or(128) * 1024 * 1024;
    }
    if let Some(v) = m.strip_suffix('K') {
        return v.parse::<usize>().unwrap_or(128) * 1024;
    }
    m.parse::<usize>().unwrap_or(128 * 1024 * 1024)
}

pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'U' => b'A',
            b'N' => b'N',
            x => x,
        })
        .collect()
}

pub fn validate_input_mode(
    command: &str,
    single_input: Option<&str>,
    paired_in1: Option<&str>,
    paired_in2: Option<&str>,
) -> Result<()> {
    if paired_in1.is_some() ^ paired_in2.is_some() {
        bail!(
            "{command}: paired-end mode requires both read files (use -i/-1 and -I/-2 together). For single-end input, pass a single positional INPUT without -i/-1."
        );
    }
    if single_input.is_some() && paired_in1.is_some() {
        bail!(
            "{command}: do not mix positional single-end INPUT with paired-end flags (-i/-1 and -I/-2)"
        );
    }
    Ok(())
}
