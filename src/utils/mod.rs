use anyhow::{bail, Result};

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

pub fn complement_in_place(seq: &mut [u8]) {
    for b in seq {
        *b = complement_base(*b);
    }
}

pub fn reverse_complement_in_place(seq: &mut [u8]) {
    seq.reverse();
    complement_in_place(seq);
}

fn complement_base(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'U' => b'A',
        b'N' => b'N',
        _ => base,
    }
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
