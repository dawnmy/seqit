use std::{
    fs::File,
    io::{BufRead, BufReader, Write},
};

use anyhow::{bail, Result};
use regex::{Regex, RegexBuilder};

use crate::cli::LocateArgs;
use crate::io;

pub fn run(args: LocateArgs) -> Result<()> {
    let pats = load_patterns(&args)?;
    let regexes = if args.regex {
        Some(
            pats.iter()
                .map(|p| {
                    RegexBuilder::new(p)
                        .case_insensitive(args.ignore_case)
                        .build()
                })
                .collect::<Result<Vec<Regex>, _>>()?,
        )
    } else {
        None
    };
    let normalized_pats = if args.ignore_case && !args.regex {
        Some(
            pats.iter()
                .map(|p| p.to_ascii_lowercase())
                .collect::<Vec<_>>(),
        )
    } else {
        None
    };
    let in_path = args.io.input.as_deref();
    let (fmt, br) = io::open_seq_reader(in_path, &args.io.format, &args.io.compression)?;
    let out = io::open_writer(&args.io.output)?;
    let out = io::wrap_compress(out, &args.io.output, &args.io.compression)?;
    let mut out = io::buffered_writer(out);
    if !matches!(
        fmt,
        crate::formats::SeqFormat::Fasta | crate::formats::SeqFormat::Fastq
    ) {
        bail!("locate currently supports FASTA/FASTQ input");
    }
    io::for_each_record_from_reader(br, fmt, |rec| {
        locate_in_record(
            rec.id_str()?,
            &rec.seq,
            &pats,
            regexes.as_deref(),
            normalized_pats.as_deref(),
            &args,
            &mut out,
        )
    })?;
    out.flush()?;
    Ok(())
}

fn locate_in_record(
    id: &str,
    seq: &[u8],
    pats: &[String],
    regexes: Option<&[Regex]>,
    normalized_pats: Option<&[String]>,
    args: &LocateArgs,
    out: &mut impl Write,
) -> Result<()> {
    let text = String::from_utf8_lossy(seq);
    let hay_lower = if args.ignore_case && !args.regex {
        Some(text.to_ascii_lowercase())
    } else {
        None
    };
    let hay = hay_lower.as_deref().unwrap_or(text.as_ref());
    for (idx, p) in pats.iter().enumerate() {
        if let Some(regexes) = regexes {
            let re = &regexes[idx];
            for m in re.find_iter(text.as_ref()) {
                emit(out, args.bed, id, m.start(), m.end(), m.as_str())?;
                if !args.all {
                    break;
                }
            }
            continue;
        }
        let query = if let Some(normalized_pats) = normalized_pats {
            normalized_pats[idx].as_str()
        } else {
            p.as_str()
        };
        let mut offset = 0usize;
        while let Some(i) = hay[offset..].find(query) {
            let s = offset + i;
            let e = s + query.len();
            emit(out, args.bed, id, s, e, &text[s..e])?;
            if !args.all {
                break;
            }
            offset = s + 1;
        }
    }
    Ok(())
}

fn emit(out: &mut impl Write, bed: bool, id: &str, s: usize, e: usize, m: &str) -> Result<()> {
    if bed {
        writeln!(out, "{id}\t{s}\t{e}\t{m}")?;
    } else {
        writeln!(out, "{id}\t{}\t{}\t{m}", s + 1, e)?;
    }
    Ok(())
}

fn load_patterns(args: &LocateArgs) -> Result<Vec<String>> {
    if let Some(p) = &args.pattern {
        let normalized = p.trim();
        if normalized.is_empty() {
            bail!("--pattern cannot be empty");
        }
        return Ok(vec![normalized.to_string()]);
    }
    if let Some(file) = &args.pattern_file {
        let f = File::open(file)?;
        let br = BufReader::with_capacity(io::BUFFER_SIZE, f);
        let patterns = br
            .lines()
            .map(|line| line.map(|l| l.trim().to_string()))
            .collect::<std::result::Result<Vec<_>, _>>()?
            .into_iter()
            .filter(|l| !l.is_empty())
            .collect::<Vec<_>>();
        if patterns.is_empty() {
            bail!("pattern file '{}' contains no non-empty patterns", file);
        }
        return Ok(patterns);
    }
    bail!("one of --pattern or --pattern-file is required")
}
