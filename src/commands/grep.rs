use std::{
    collections::HashSet,
    fs::File,
    io::{BufRead, BufReader},
};

use anyhow::{bail, Context, Result};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use regex::{Regex, RegexBuilder};

use crate::{
    cli::{GrepArgs, SearchBy},
    io, utils,
};

pub fn run(args: GrepArgs) -> Result<()> {
    let patterns = load_patterns(&args)?;
    if !args.io.quiet {
        if let Some(path) = &args.pattern_file {
            eprintln!(
                "[INFO] {} patterns loaded from file '{}'",
                patterns.len(),
                path
            );
        } else {
            eprintln!("[INFO] {} pattern(s) loaded", patterns.len());
        }
        if args.count {
            eprintln!("[INFO] grep count mode");
        } else if args.only_names {
            eprintln!("[INFO] grep names mode");
        } else {
            eprintln!("[INFO] grep record output mode");
        }
    }
    let matcher = Matcher::new(&patterns, &args.by, &args)?;
    utils::validate_input_mode(
        "grep",
        args.io.input.as_deref(),
        args.in1.as_deref(),
        args.in2.as_deref(),
    )?;

    if let (Some(in1), Some(in2)) = (args.in1.as_deref(), args.in2.as_deref()) {
        let out2 = args
            .output2
            .as_deref()
            .context("paired grep requires --output2")?;
        let mut w1 = None;
        let mut w2 = None;
        if !args.count && !args.only_names {
            let o1 = io::open_writer(&args.io.output)?;
            let o1 = io::wrap_compress_for_streams(o1, &args.io.output, &args.io.compression, 2)?;
            let o2 = io::open_writer(out2)?;
            let o2 = io::wrap_compress_for_streams(o2, out2, &args.io.compression, 2)?;
            w1 = Some(io::buffered_writer(o1));
            w2 = Some(io::buffered_writer(o2));
        }
        let mut count = 0usize;
        let mut invalid_preview = Vec::new();
        let mut invalid_count = 0usize;
        let progress = make_progress_bar(args.progress);
        io::for_each_fastq_pair(in1, in2, &args.io.compression, |pair| {
            match pair {
                io::FastqPair::Both(a, b) => {
                    let a_id = a.id_str()?;
                    let b_id = b.id_str()?;
                    let a_desc = a.desc_str()?;
                    let b_desc = b.desc_str()?;
                    if pair_key(a_id) != pair_key(b_id) {
                        invalid_count += 1;
                        if invalid_preview.len() < 10 {
                            invalid_preview.push(format!("R1={a_id} R2={b_id}"));
                        }
                        return Ok(true);
                    }
                    let mut keep = if args.pair_mode == "both" {
                        is_match_fields(a_id, a_desc, &a.seq, a.qual, &matcher, &args.by)
                            && is_match_fields(b_id, b_desc, &b.seq, b.qual, &matcher, &args.by)
                    } else {
                        is_match_fields(a_id, a_desc, &a.seq, a.qual, &matcher, &args.by)
                            || is_match_fields(b_id, b_desc, &b.seq, b.qual, &matcher, &args.by)
                    };
                    if args.invert {
                        keep = !keep;
                    }
                    if keep {
                        count += 1;
                        if args.only_names {
                            println!("{a_id}");
                        } else if let (Some(w1), Some(w2)) = (&mut w1, &mut w2) {
                            io::write_record_ref(w1, crate::formats::SeqFormat::Fastq, &a)?;
                            io::write_record_ref(w2, crate::formats::SeqFormat::Fastq, &b)?;
                        }
                    }
                }
                io::FastqPair::Left(a) => {
                    invalid_count += 1;
                    if invalid_preview.len() < 10 {
                        invalid_preview.push(format!("unmatched in R1: {}", a.id_str()?));
                    }
                }
                io::FastqPair::Right(b) => {
                    invalid_count += 1;
                    if invalid_preview.len() < 10 {
                        invalid_preview.push(format!("unmatched in R2: {}", b.id_str()?));
                    }
                }
            }
            if let Some(pb) = &progress {
                pb.inc(1);
            }
            Ok(true)
        })?;
        report_invalid_pairs(invalid_count, invalid_preview, args.allow_unpaired)?;
        if let Some(pb) = &progress {
            pb.finish_and_clear();
        }
        if args.count {
            println!("{count}");
        }
        return Ok(());
    }

    let in_path = args.io.input.as_deref();
    let (fmt, br) = io::open_seq_reader(in_path, &args.io.format, &args.io.compression)?;
    let progress = make_progress_bar(args.progress);
    let mut count = 0usize;
    let mut writer = None;
    if !args.count && !args.only_names {
        let out = io::open_writer(&args.io.output)?;
        let out = io::wrap_compress(out, &args.io.output, &args.io.compression)?;
        writer = Some(io::buffered_writer(out));
    }
    io::for_each_record_from_reader(br, fmt, |rec| {
        let id = rec.id_str()?;
        let desc = rec.desc_str()?;
        let mut keep = is_match_fields(id, desc, &rec.seq, rec.qual, &matcher, &args.by);
        if args.invert {
            keep = !keep;
        }
        if keep {
            count += 1;
            if args.only_names {
                println!("{id}");
            } else if let Some(w) = &mut writer {
                io::write_record_ref(w, fmt, &rec)?;
            }
        }
        if let Some(pb) = &progress {
            pb.inc(1);
        }
        Ok(())
    })?;
    if let Some(pb) = &progress {
        pb.finish_and_clear();
    }
    if args.count {
        println!("{count}");
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

fn make_progress_bar(enabled: bool) -> Option<ProgressBar> {
    if !enabled {
        return None;
    }
    let pb = ProgressBar::with_draw_target(None, ProgressDrawTarget::stderr());
    pb.set_style(
        ProgressStyle::with_template("{spinner:.green} processed {pos} records")
            .expect("valid progress template"),
    );
    Some(pb)
}

fn is_match_fields(
    id: &str,
    desc: Option<&str>,
    seq: &[u8],
    qual: Option<&[u8]>,
    matcher: &Matcher,
    by: &SearchBy,
) -> bool {
    match by {
        SearchBy::Id => matcher.is_match(id),
        SearchBy::Name => {
            if let Some(desc) = desc {
                matcher.is_match(&format!("{id} {desc}"))
            } else {
                matcher.is_match(id)
            }
        }
        SearchBy::Seq => {
            let seq = String::from_utf8_lossy(seq);
            matcher.is_match(&seq)
        }
        SearchBy::Qual => qual
            .map(|q| {
                let qual = String::from_utf8_lossy(q);
                matcher.is_match(&qual)
            })
            .unwrap_or(false),
    }
}

fn pair_key(id: &str) -> &str {
    id.strip_suffix("/1")
        .or_else(|| id.strip_suffix("/2"))
        .unwrap_or(id)
}

fn report_invalid_pairs(
    invalid_count: usize,
    invalid_preview: Vec<String>,
    allow_unpaired: bool,
) -> Result<()> {
    if invalid_count == 0 {
        return Ok(());
    }
    let preview = invalid_preview.join("; ");
    if !allow_unpaired {
        bail!(
            "paired inputs contain {invalid_count} invalid/unpaired records; examples: {preview}. rerun with --allow-unpaired to continue"
        );
    }
    eprintln!(
        "warning: paired inputs contain {invalid_count} invalid/unpaired records; examples: {preview}. unpaired records will be skipped"
    );
    Ok(())
}
