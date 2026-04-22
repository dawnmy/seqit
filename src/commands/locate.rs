use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use anyhow::{bail, Result};
use regex::RegexBuilder;

use crate::{cli::LocateArgs, formats::SeqFormat, io};

pub fn run(args: LocateArgs) -> Result<()> {
    let pats = load_patterns(&args)?;
    let in_path = args.io.input.as_deref();
    let fmt = SeqFormat::from_arg(&args.io.format).unwrap_or(SeqFormat::detect(in_path)?);
    let recs = io::read_records(in_path, fmt, &args.io.compression)?;

    for r in recs {
        let text = String::from_utf8_lossy(&r.seq).to_string();
        for p in &pats {
            if args.regex {
                let re = RegexBuilder::new(p)
                    .case_insensitive(args.ignore_case)
                    .build()?;
                for m in re.find_iter(&text) {
                    emit(args.bed, &r.id, m.start(), m.end(), m.as_str());
                    if !args.all {
                        break;
                    }
                }
            } else {
                let query = if args.ignore_case {
                    p.to_ascii_lowercase()
                } else {
                    p.clone()
                };
                let hay = if args.ignore_case {
                    text.to_ascii_lowercase()
                } else {
                    text.clone()
                };
                let mut offset = 0usize;
                while let Some(i) = hay[offset..].find(&query) {
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
        return Ok(vec![p.clone()]);
    }
    if let Some(file) = &args.pattern_file {
        let f = File::open(file)?;
        let br = BufReader::new(f);
        return Ok(br.lines().collect::<std::result::Result<Vec<_>, _>>()?);
    }
    bail!("one of --pattern or --pattern-file is required")
}
