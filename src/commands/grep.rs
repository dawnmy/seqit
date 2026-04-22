use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use anyhow::{bail, Context, Result};
use regex::RegexBuilder;

use crate::{
    cli::{GrepArgs, SearchBy},
    formats::{SeqFormat, SeqRecord},
    io, pairs,
};

pub fn run(args: GrepArgs) -> Result<()> {
    let patterns = load_patterns(&args)?;

    if let (Some(in1), Some(in2)) = (args.in1.as_deref(), args.in2.as_deref()) {
        let r1 = io::read_records(Some(in1), SeqFormat::Fastq, &args.io.compression)?;
        let r2 = io::read_records(Some(in2), SeqFormat::Fastq, &args.io.compression)?;
        pairs::validate_pair_counts(&r1, &r2)?;
        let mut o1 = Vec::new();
        let mut o2 = Vec::new();
        let mut count = 0usize;
        for (a, b) in r1.into_iter().zip(r2.into_iter()) {
            let m1 = is_match(&a, &patterns, &args)?;
            let m2 = is_match(&b, &patterns, &args)?;
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
    let mut out = Vec::new();
    let mut count = 0usize;
    for r in recs {
        let mut keep = is_match(&r, &patterns, &args)?;
        if args.invert {
            keep = !keep;
        }
        if keep {
            count += 1;
            if !args.count {
                out.push(r);
            }
        }
    }
    if args.count {
        println!("{count}");
    } else if args.only_names {
        for r in out {
            println!("{}", r.id);
        }
    } else {
        io::write_records(&args.io.output, fmt, &args.io.compression, &out)?;
    }
    Ok(())
}

fn load_patterns(args: &GrepArgs) -> Result<Vec<String>> {
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

fn is_match(rec: &SeqRecord, patterns: &[String], args: &GrepArgs) -> Result<bool> {
    let target = match args.by {
        SearchBy::Id | SearchBy::Name => rec.name(),
        SearchBy::Seq => String::from_utf8_lossy(&rec.seq).to_string(),
        SearchBy::Qual => rec
            .qual
            .as_ref()
            .map(|q| String::from_utf8_lossy(q).to_string())
            .unwrap_or_default(),
    };
    if args.regex {
        for p in patterns {
            let re = RegexBuilder::new(p)
                .case_insensitive(args.ignore_case)
                .build()?;
            if re.is_match(&target) {
                return Ok(true);
            }
        }
        Ok(false)
    } else {
        if args.ignore_case {
            let t = target.to_ascii_lowercase();
            Ok(patterns.iter().any(|p| t.contains(&p.to_ascii_lowercase())))
        } else {
            Ok(patterns.iter().any(|p| target.contains(p)))
        }
    }
}
