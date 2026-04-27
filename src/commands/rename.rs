use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, Write},
};

use anyhow::{bail, Context, Result};
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
    let mut out = io::buffered_writer(w);
    let mut idx = 0usize;
    if !matches!(fmt, SeqFormat::Fasta | SeqFormat::Fastq) {
        bail!("rename currently supports FASTA/FASTQ input");
    }
    io::for_each_record_from_reader(br, fmt, |rec| {
        let id = rename_id(args, mode, idx, rec.id_str()?, mapping, regex, replacement)?;
        match fmt {
            SeqFormat::Fasta => {
                out.write_all(b">")?;
                out.write_all(id.as_bytes())?;
                if let Some(desc) = rec.desc {
                    out.write_all(b" ")?;
                    out.write_all(desc)?;
                }
                out.write_all(b"\n")?;
                out.write_all(&rec.seq)?;
                out.write_all(b"\n")?;
            }
            SeqFormat::Fastq => {
                let Some(qual) = rec.qual else {
                    bail!(
                        "record '{}' missing qualities for FASTQ output",
                        rec.id_str()?
                    );
                };
                out.write_all(b"@")?;
                out.write_all(id.as_bytes())?;
                if let Some(desc) = rec.desc {
                    out.write_all(b" ")?;
                    out.write_all(desc)?;
                }
                out.write_all(b"\n")?;
                out.write_all(&rec.seq)?;
                out.write_all(b"\n+\n")?;
                out.write_all(qual)?;
                out.write_all(b"\n")?;
            }
            _ => unreachable!("format checked before streaming"),
        }
        idx += 1;
        Ok(())
    })?;
    out.flush()?;
    Ok(())
}

fn run_paired_streaming(
    args: &RenameArgs,
    in1: &str,
    in2: &str,
    out2: &str,
    plan: RenamePlan<'_>,
) -> Result<()> {
    let w1 = io::open_writer(&args.io.output)?;
    let w1 = io::wrap_compress_for_streams(w1, &args.io.output, &args.io.compression, 2)?;
    let w2 = io::open_writer(out2)?;
    let w2 = io::wrap_compress_for_streams(w2, out2, &args.io.compression, 2)?;
    let mut w1 = io::buffered_writer(w1);
    let mut w2 = io::buffered_writer(w2);

    let mut idx = 0usize;
    let mut invalid_preview = Vec::new();
    let mut invalid_count = 0usize;
    io::for_each_fastq_pair(in1, in2, &args.io.compression, |pair| {
        match pair {
            io::FastqPair::Both(a, b) => {
                let a_id = a.id_str()?;
                let b_id = b.id_str()?;
                if pair_key(a_id) != pair_key(b_id) {
                    invalid_count += 1;
                    if invalid_preview.len() < 10 {
                        invalid_preview.push(format!("R1={a_id} R2={b_id}"));
                    }
                    return Ok(true);
                }
                let core = rename_id(
                    args,
                    plan.mode,
                    idx,
                    a_id,
                    plan.mapping,
                    plan.regex,
                    plan.replacement,
                )?;
                let (id1, id2) = if args.keep_pair_suffix {
                    (format!("{core}/1"), format!("{core}/2"))
                } else {
                    (core.clone(), core)
                };
                write_renamed_fastq(&mut w1, &id1, a.desc, &a.seq, a.qual)?;
                write_renamed_fastq(&mut w2, &id2, b.desc, &b.seq, b.qual)?;
                idx += 1;
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
        Ok(true)
    })?;
    report_invalid_pairs(invalid_count, invalid_preview, args.allow_unpaired)
}

fn write_renamed_fastq(
    out: &mut impl Write,
    id: &str,
    desc: Option<&[u8]>,
    seq: &[u8],
    qual: Option<&[u8]>,
) -> Result<()> {
    let Some(qual) = qual else {
        bail!("record '{id}' missing qualities for FASTQ output");
    };
    out.write_all(b"@")?;
    out.write_all(id.as_bytes())?;
    if let Some(desc) = desc {
        out.write_all(b" ")?;
        out.write_all(desc)?;
    }
    out.write_all(b"\n")?;
    out.write_all(seq)?;
    out.write_all(b"\n+\n")?;
    out.write_all(qual)?;
    out.write_all(b"\n")?;
    Ok(())
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
