use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use anyhow::{bail, Result};
use regex::{Regex, RegexBuilder};

use crate::{cli::LocateArgs, formats::SeqFormat, io};

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
    let fmt = SeqFormat::from_arg_or_detect(&args.io.format, in_path)?;
    let recs = io::read_records(in_path, fmt, &args.io.compression)?;

    for r in recs {
        let text = String::from_utf8_lossy(&r.seq).to_string();
        for (idx, p) in pats.iter().enumerate() {
            if let Some(regexes) = &regexes {
                let re = &regexes[idx];
                for m in re.find_iter(&text) {
                    emit(args.bed, &r.id, m.start(), m.end(), m.as_str());
                    if !args.all {
                        break;
                    }
                }
            } else {
                let query = if let Some(normalized_pats) = &normalized_pats {
                    normalized_pats[idx].as_str()
                } else {
                    p.as_str()
                };
                let hay = if args.ignore_case {
                    text.to_ascii_lowercase()
                } else {
                    text.clone()
                };
                let mut offset = 0usize;
                while let Some(i) = hay[offset..].find(query) {
                    let s = offset + i;
                    let e = s + query.len();
                    emit(args.bed, &r.id, s, e, &text[s..e]);
                    if !args.all {
                        break;
                    }
                    offset = s + 1;
                }
            }
        }
    }
    Ok(())
}

fn emit(bed: bool, id: &str, s: usize, e: usize, m: &str) {
    if bed {
        println!("{id}\t{s}\t{e}\t{m}");
    } else {
        println!("{id}\t{}\t{}\t{m}", s + 1, e);
    }
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
        let br = BufReader::new(f);
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
