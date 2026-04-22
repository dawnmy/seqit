use std::{
    collections::HashSet,
    fs::File,
    io::{BufRead, BufReader},
};

use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use regex::{Regex, RegexBuilder};

use crate::{
    cli::{GrepArgs, SearchBy},
    formats::{SeqFormat, SeqRecord},
    io, pairs,
};

pub fn run(args: GrepArgs) -> Result<()> {
    let patterns = load_patterns(&args)?;
    let matcher = Matcher::new(&patterns, &args.by, &args)?;

    if let (Some(in1), Some(in2)) = (args.in1.as_deref(), args.in2.as_deref()) {
        let r1 = io::read_records(Some(in1), SeqFormat::Fastq, &args.io.compression)?;
        let r2 = io::read_records(Some(in2), SeqFormat::Fastq, &args.io.compression)?;
        let (r1, r2) = pairs::prepare_paired_records(r1, r2, args.allow_unpaired)?;
        let mut o1 = Vec::new();
        let mut o2 = Vec::new();
        let mut count = 0usize;
        let total = r1.len().min(r2.len());
        for (idx, (a, b)) in r1.into_iter().zip(r2.into_iter()).enumerate() {
            let m1 = is_match(&a, &matcher, &args.by);
            let m2 = is_match(&b, &matcher, &args.by);
            let mut keep = if args.pair_mode == "both" {
                m1 && m2
            } else {
                m1 || m2
            };
            if args.invert {
                keep = !keep;
            }
            if keep {
                count += 1;
                if !args.count {
                    o1.push(a);
                    o2.push(b);
                }
            }
            report_progress(args.progress, idx + 1, total);
        }
        if args.count {
            println!("{count}");
            return Ok(());
        }
        let out2 = args
            .output2
            .as_deref()
            .context("paired grep requires --output2")?;
        if args.only_names {
            for r in &o1 {
                println!("{}", r.id);
            }
            return Ok(());
        }
        io::write_records(&args.io.output, SeqFormat::Fastq, &args.io.compression, &o1)?;
        io::write_records(out2, SeqFormat::Fastq, &args.io.compression, &o2)?;
        return Ok(());
    }

    let in_path = args.io.input.as_deref();
    let fmt = SeqFormat::from_arg(&args.io.format).unwrap_or(SeqFormat::detect(in_path)?);
    let recs = io::read_records(in_path, fmt, &args.io.compression)?;
    let selected: Vec<&SeqRecord> = if args.progress {
        let mut out = Vec::new();
        let total = recs.len();
        for (idx, rec) in recs.iter().enumerate() {
            let mut keep = is_match(rec, &matcher, &args.by);
            if args.invert {
                keep = !keep;
            }
            if keep {
                out.push(rec);
            }
            report_progress(true, idx + 1, total);
        }
        out
    } else {
        recs.par_iter()
            .filter(|r| {
                let mut keep = is_match(r, &matcher, &args.by);
                if args.invert {
                    keep = !keep;
                }
                keep
            })
            .collect()
    };

    if args.count {
        println!("{}", selected.len());
    } else if args.only_names {
        for r in selected {
            println!("{}", r.id);
        }
    } else {
        let out: Vec<SeqRecord> = selected.into_iter().cloned().collect();
        io::write_records(&args.io.output, fmt, &args.io.compression, &out)?;
    }
    Ok(())
}

fn load_patterns(args: &GrepArgs) -> Result<Vec<String>> {
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

enum Matcher {
    Regex(Vec<Regex>),
    ExactCaseSensitive(HashSet<String>),
    ExactCaseInsensitive(HashSet<String>),
    CaseSensitive(Vec<String>),
    CaseInsensitive(Vec<String>),
}

impl Matcher {
    fn new(patterns: &[String], by: &SearchBy, args: &GrepArgs) -> Result<Self> {
        if args.regex {
            let regexes = patterns
                .iter()
                .map(|p| {
                    RegexBuilder::new(p)
                        .case_insensitive(args.ignore_case)
                        .build()
                })
                .collect::<Result<Vec<_>, _>>()?;
            Ok(Self::Regex(regexes))
        } else if matches!(by, SearchBy::Id) && args.ignore_case {
            Ok(Self::ExactCaseInsensitive(
                patterns.iter().map(|p| p.to_ascii_lowercase()).collect(),
            ))
        } else if matches!(by, SearchBy::Id) {
            Ok(Self::ExactCaseSensitive(patterns.iter().cloned().collect()))
        } else if args.ignore_case {
            Ok(Self::CaseInsensitive(
                patterns.iter().map(|p| p.to_ascii_lowercase()).collect(),
            ))
        } else {
            Ok(Self::CaseSensitive(patterns.to_vec()))
        }
    }

    fn is_match(&self, target: &str) -> bool {
        match self {
            Self::Regex(regexes) => regexes.iter().any(|re| re.is_match(target)),
            Self::ExactCaseSensitive(patterns) => patterns.contains(target),
            Self::ExactCaseInsensitive(patterns) => {
                let lowered = target.to_ascii_lowercase();
                patterns.contains(&lowered)
            }
            Self::CaseSensitive(patterns) => patterns.iter().any(|p| target.contains(p)),
            Self::CaseInsensitive(patterns) => {
                let lowered = target.to_ascii_lowercase();
                patterns.iter().any(|p| lowered.contains(p))
            }
        }
    }
}

fn report_progress(enabled: bool, processed: usize, total: usize) {
    if !enabled || total == 0 {
        return;
    }
    const STEP: usize = 100_000;
    if processed == total || processed % STEP == 0 {
        eprintln!(
            "grep progress: {processed}/{total} ({:.1}%)",
            processed as f64 * 100.0 / total as f64
        );
    }
}

fn is_match(rec: &SeqRecord, matcher: &Matcher, by: &SearchBy) -> bool {
    match by {
        SearchBy::Id => matcher.is_match(&rec.id),
        SearchBy::Name => matcher.is_match(&rec.name()),
        SearchBy::Seq => {
            let seq = String::from_utf8_lossy(&rec.seq);
            matcher.is_match(&seq)
        }
        SearchBy::Qual => rec
            .qual
            .as_ref()
            .map(|q| {
                let qual = String::from_utf8_lossy(q);
                matcher.is_match(&qual)
            })
            .unwrap_or(false),
    }
}
