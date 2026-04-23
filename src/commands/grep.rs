use std::{
    collections::HashSet,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
};

use anyhow::{bail, Context, Result};
use bio::io::{fasta, fastq};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use regex::{Regex, RegexBuilder};

use crate::{
    cli::{GrepArgs, SearchBy},
    formats::SeqFormat,
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
        let r1 = io::open_reader(Some(in1))?;
        let r1 = io::wrap_decompress(r1, Some(in1), &args.io.compression)?;
        let r2 = io::open_reader(Some(in2))?;
        let r2 = io::wrap_decompress(r2, Some(in2), &args.io.compression)?;
        let fq1 = fastq::Reader::new(BufReader::new(r1));
        let fq2 = fastq::Reader::new(BufReader::new(r2));
        let mut it1 = fq1.records();
        let mut it2 = fq2.records();

        let mut w1 = None;
        let mut w2 = None;
        if !args.count && !args.only_names {
            let o1 = io::open_writer(&args.io.output)?;
            let o1 = io::wrap_compress(o1, &args.io.output, &args.io.compression)?;
            let o2 = io::open_writer(out2)?;
            let o2 = io::wrap_compress(o2, out2, &args.io.compression)?;
            w1 = Some(fastq::Writer::new(BufWriter::new(o1)));
            w2 = Some(fastq::Writer::new(BufWriter::new(o2)));
        }
        let mut count = 0usize;
        let mut invalid_preview = Vec::new();
        let mut invalid_count = 0usize;
        let progress = make_progress_bar(args.progress);
        loop {
            let a = it1.next().transpose()?;
            let b = it2.next().transpose()?;
            match (a, b) {
                (Some(a), Some(b)) => {
                    if pair_key(a.id()) != pair_key(b.id()) {
                        invalid_count += 1;
                        if invalid_preview.len() < 10 {
                            invalid_preview.push(format!("R1={} R2={}", a.id(), b.id()));
                        }
                        continue;
                    }
                    let mut keep = if args.pair_mode == "both" {
                        is_match_fields(
                            a.id(),
                            a.desc(),
                            a.seq(),
                            Some(a.qual()),
                            &matcher,
                            &args.by,
                        ) && is_match_fields(
                            b.id(),
                            b.desc(),
                            b.seq(),
                            Some(b.qual()),
                            &matcher,
                            &args.by,
                        )
                    } else {
                        is_match_fields(
                            a.id(),
                            a.desc(),
                            a.seq(),
                            Some(a.qual()),
                            &matcher,
                            &args.by,
                        ) || is_match_fields(
                            b.id(),
                            b.desc(),
                            b.seq(),
                            Some(b.qual()),
                            &matcher,
                            &args.by,
                        )
                    };
                    if args.invert {
                        keep = !keep;
                    }
                    if keep {
                        count += 1;
                        if args.only_names {
                            println!("{}", a.id());
                        } else if let (Some(w1), Some(w2)) = (&mut w1, &mut w2) {
                            w1.write(a.id(), a.desc(), a.seq(), a.qual())?;
                            w2.write(b.id(), b.desc(), b.seq(), b.qual())?;
                        }
                    }
                }
                (Some(a), None) => {
                    invalid_count += 1;
                    if invalid_preview.len() < 10 {
                        invalid_preview.push(format!("unmatched in R1: {}", a.id()));
                    }
                }
                (None, Some(b)) => {
                    invalid_count += 1;
                    if invalid_preview.len() < 10 {
                        invalid_preview.push(format!("unmatched in R2: {}", b.id()));
                    }
                }
                (None, None) => break,
            }
            if let Some(pb) = &progress {
                pb.inc(1);
            }
        }
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
        writer = Some(BufWriter::new(out));
    }
    match fmt {
        SeqFormat::Fasta => {
            let fa = fasta::Reader::new(br);
            for rec in fa.records() {
                let rec = rec?;
                let mut keep =
                    is_match_fields(rec.id(), rec.desc(), rec.seq(), None, &matcher, &args.by);
                if args.invert {
                    keep = !keep;
                }
                if keep {
                    count += 1;
                    if args.only_names {
                        println!("{}", rec.id());
                    } else if let Some(w) = &mut writer {
                        write_fasta_record(w, rec.id(), rec.desc(), rec.seq())?;
                    }
                }
                if let Some(pb) = &progress {
                    pb.inc(1);
                }
            }
        }
        SeqFormat::Fastq => {
            let fq = fastq::Reader::new(br);
            for rec in fq.records() {
                let rec = rec?;
                let mut keep = is_match_fields(
                    rec.id(),
                    rec.desc(),
                    rec.seq(),
                    Some(rec.qual()),
                    &matcher,
                    &args.by,
                );
                if args.invert {
                    keep = !keep;
                }
                if keep {
                    count += 1;
                    if args.only_names {
                        println!("{}", rec.id());
                    } else if let Some(w) = &mut writer {
                        write_fastq_record(w, rec.id(), rec.desc(), rec.seq(), rec.qual())?;
                    }
                }
                if let Some(pb) = &progress {
                    pb.inc(1);
                }
            }
        }
        _ => bail!("grep currently supports FASTA/FASTQ input"),
    }
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

fn write_fasta_record(w: &mut impl Write, id: &str, desc: Option<&str>, seq: &[u8]) -> Result<()> {
    if let Some(desc) = desc {
        writeln!(w, ">{id} {desc}")?;
    } else {
        writeln!(w, ">{id}")?;
    }
    writeln!(w, "{}", String::from_utf8_lossy(seq))?;
    Ok(())
}

fn write_fastq_record(
    w: &mut impl Write,
    id: &str,
    desc: Option<&str>,
    seq: &[u8],
    qual: &[u8],
) -> Result<()> {
    if let Some(desc) = desc {
        writeln!(w, "@{id} {desc}")?;
    } else {
        writeln!(w, "@{id}")?;
    }
    writeln!(w, "{}", String::from_utf8_lossy(seq))?;
    writeln!(w, "+")?;
    writeln!(w, "{}", String::from_utf8_lossy(qual))?;
    Ok(())
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
