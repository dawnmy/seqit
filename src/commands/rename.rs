use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
};

use anyhow::{bail, Context, Result};
use bio::io::{fasta, fastq};
use regex::Regex;

use crate::{
    cli::{RenameArgs, RenameMode},
    formats::SeqFormat,
    io, utils,
};

pub fn run(args: RenameArgs) -> Result<()> {
    let mode = resolve_mode(&args)?;
    let mapping = if mode == RenameMode::Map {
        Some(load_mapping(
            args.map_file.as_deref().context("missing --map-file")?,
        )?)
    } else {
        None
    };
    let regex = if mode == RenameMode::Regex {
        let pattern = args
            .match_regex
            .as_deref()
            .context("missing --match-regex")?;
        Some(Regex::new(pattern).context("invalid --match-regex pattern")?)
    } else {
        None
    };
    let replacement = if mode == RenameMode::Regex {
        Some(args.replace.as_deref().context("missing --replace")?)
    } else {
        None
    };

    let paired_in1 = args.input1.as_deref().or(args.in1.as_deref());
    let paired_in2 = args.input2.as_deref().or(args.in2.as_deref());
    utils::validate_input_mode("rename", args.io.input.as_deref(), paired_in1, paired_in2)?;
    if let (Some(in1), Some(in2)) = (paired_in1, paired_in2) {
        let out2 = args
            .output2
            .as_deref()
            .context("paired rename requires --output2")?;
        return run_paired_streaming(
            &args,
            in1,
            in2,
            out2,
            RenamePlan {
                mode,
                mapping: mapping.as_ref(),
                regex: regex.as_ref(),
                replacement,
            },
        );
    }

    let in_path = args.io.input.as_deref();
    run_single_streaming(
        &args,
        in_path,
        mode,
        mapping.as_ref(),
        regex.as_ref(),
        replacement,
    )
}

fn resolve_mode(args: &RenameArgs) -> Result<RenameMode> {
    let has_map = args.map_file.is_some();
    let has_regex = args.match_regex.is_some() || args.replace.is_some();
    let has_generate_opts =
        args.template.is_some() || args.prefix != "seq" || args.start != 1 || args.width != 6;

    if has_map && has_regex {
        bail!("--map-file cannot be combined with --match-regex/--replace");
    }
    match args.mode {
        RenameMode::Auto => {
            if has_map {
                Ok(RenameMode::Map)
            } else if has_regex {
                Ok(RenameMode::Regex)
            } else {
                Ok(RenameMode::Generate)
            }
        }
        RenameMode::Generate => {
            if has_map || has_regex {
                bail!("--mode generate cannot be combined with mapping or regex options");
            }
            Ok(RenameMode::Generate)
        }
        RenameMode::Map => {
            if args.map_file.is_none() {
                bail!("--mode map requires --map-file");
            }
            if has_regex || has_generate_opts {
                bail!("--mode map cannot be combined with generation or regex options");
            }
            Ok(RenameMode::Map)
        }
        RenameMode::Regex => {
            if args.match_regex.is_none() || args.replace.is_none() {
                bail!("--mode regex requires both --match-regex and --replace");
            }
            if has_map || has_generate_opts {
                bail!("--mode regex cannot be combined with generation or mapping options");
            }
            Ok(RenameMode::Regex)
        }
    }
}

fn load_mapping(path: &str) -> Result<HashMap<String, String>> {
    let mut map = HashMap::new();
    let file = File::open(path).with_context(|| format!("failed to read mapping file '{path}'"))?;
    let br = BufReader::with_capacity(io::BUFFER_SIZE, file);
    for (i, raw_line) in br.lines().enumerate() {
        let line_no = i + 1;
        let line_raw = raw_line?;
        let line = line_raw.trim_end();
        if line.is_empty() {
            continue;
        }
        let mut cols = line.splitn(3, '\t');
        let old = cols
            .next()
            .filter(|c| !c.is_empty())
            .with_context(|| format!("invalid mapping line {line_no}: missing old ID"))?;
        let new = cols
            .next()
            .filter(|c| !c.is_empty())
            .with_context(|| format!("invalid mapping line {line_no}: missing new ID"))?;
        if cols.next().is_some() {
            bail!("invalid mapping line {line_no}: expected exactly two tab-separated columns");
        }
        map.insert(old.to_string(), new.to_string());
    }
    Ok(map)
}

fn rename_id(
    args: &RenameArgs,
    mode: RenameMode,
    idx: usize,
    old_id: &str,
    mapping: Option<&HashMap<String, String>>,
    regex: Option<&Regex>,
    replacement: Option<&str>,
) -> Result<String> {
    match mode {
        RenameMode::Generate => Ok(format_name(args, args.start + idx)),
        RenameMode::Map => Ok(mapping
            .and_then(|m| m.get(old_id))
            .cloned()
            .unwrap_or_else(|| old_id.to_string())),
        RenameMode::Regex => Ok(regex
            .context("missing compiled regex")?
            .replace_all(old_id, replacement.context("missing replacement")?)
            .to_string()),
        RenameMode::Auto => unreachable!("mode must be resolved before renaming"),
    }
}

fn format_name(args: &RenameArgs, n: usize) -> String {
    if let Some(tpl) = &args.template {
        tpl.replace("{prefix}", &args.prefix)
            .replace("{n}", &format!("{:0width$}", n, width = args.width))
    } else {
        format!("{}{:0width$}", args.prefix, n, width = args.width)
    }
}

fn run_single_streaming(
    args: &RenameArgs,
    in_path: Option<&str>,
    mode: RenameMode,
    mapping: Option<&HashMap<String, String>>,
    regex: Option<&Regex>,
    replacement: Option<&str>,
) -> Result<()> {
    let (fmt, br) = io::open_seq_reader(in_path, &args.io.format, &args.io.compression)?;
    let w = io::open_writer(&args.io.output)?;
    let w = io::wrap_compress(w, &args.io.output, &args.io.compression)?;
    let mut idx = 0usize;
    match fmt {
        SeqFormat::Fasta => {
            let fa = fasta::Reader::new(br);
            let mut out = fasta::Writer::new(io::buffered_writer(w));
            for rec in fa.records() {
                let rec = rec?;
                let id = rename_id(args, mode, idx, rec.id(), mapping, regex, replacement)?;
                out.write(&id, rec.desc(), rec.seq())?;
                idx += 1;
            }
        }
        SeqFormat::Fastq => {
            let fq = fastq::Reader::new(br);
            let mut out = fastq::Writer::new(io::buffered_writer(w));
            for rec in fq.records() {
                let rec = rec?;
                let id = rename_id(args, mode, idx, rec.id(), mapping, regex, replacement)?;
                out.write(&id, rec.desc(), rec.seq(), rec.qual())?;
                idx += 1;
            }
        }
        _ => bail!("rename currently supports FASTA/FASTQ input"),
    }
    Ok(())
}

fn run_paired_streaming(
    args: &RenameArgs,
    in1: &str,
    in2: &str,
    out2: &str,
    plan: RenamePlan<'_>,
) -> Result<()> {
    let r1 = io::open_reader(Some(in1))?;
    let r1 = io::wrap_decompress(r1, Some(in1), &args.io.compression)?;
    let r2 = io::open_reader(Some(in2))?;
    let r2 = io::wrap_decompress(r2, Some(in2), &args.io.compression)?;
    let fq1 = fastq::Reader::new(io::buffered_reader(r1));
    let fq2 = fastq::Reader::new(io::buffered_reader(r2));
    let mut it1 = fq1.records();
    let mut it2 = fq2.records();

    let w1 = io::open_writer(&args.io.output)?;
    let w1 = io::wrap_compress_for_streams(w1, &args.io.output, &args.io.compression, 2)?;
    let w2 = io::open_writer(out2)?;
    let w2 = io::wrap_compress_for_streams(w2, out2, &args.io.compression, 2)?;
    let mut w1 = fastq::Writer::new(io::buffered_writer(w1));
    let mut w2 = fastq::Writer::new(io::buffered_writer(w2));

    let mut idx = 0usize;
    let mut invalid_preview = Vec::new();
    let mut invalid_count = 0usize;
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
                let core = rename_id(
                    args,
                    plan.mode,
                    idx,
                    a.id(),
                    plan.mapping,
                    plan.regex,
                    plan.replacement,
                )?;
                let (id1, id2) = if args.keep_pair_suffix {
                    (format!("{core}/1"), format!("{core}/2"))
                } else {
                    (core.clone(), core)
                };
                w1.write(&id1, a.desc(), a.seq(), a.qual())?;
                w2.write(&id2, b.desc(), b.seq(), b.qual())?;
                idx += 1;
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
    }
    report_invalid_pairs(invalid_count, invalid_preview, args.allow_unpaired)
}

struct RenamePlan<'a> {
    mode: RenameMode,
    mapping: Option<&'a HashMap<String, String>>,
    regex: Option<&'a Regex>,
    replacement: Option<&'a str>,
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
