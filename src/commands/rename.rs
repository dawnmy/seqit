use std::{collections::HashMap, fs};

use anyhow::{bail, Context, Result};
use regex::Regex;

use crate::{
    cli::{RenameArgs, RenameMode},
    formats::SeqFormat,
    io, pairs,
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
    if let (Some(in1), Some(in2)) = (paired_in1, paired_in2) {
        let r1 = io::read_records(Some(in1), SeqFormat::Fastq, &args.io.compression)?;
        let r2 = io::read_records(Some(in2), SeqFormat::Fastq, &args.io.compression)?;
        let (mut r1, mut r2) = pairs::prepare_paired_records(r1, r2, args.allow_unpaired)?;
        let out2 = args
            .output2
            .as_deref()
            .context("paired rename requires --output2")?;
        for (i, (a, b)) in r1.iter_mut().zip(r2.iter_mut()).enumerate() {
            let core = rename_id(&args, mode, i, &a.id, &mapping, regex.as_ref(), replacement)?;
            if args.keep_pair_suffix {
                a.id = format!("{core}/1");
                b.id = format!("{core}/2");
            } else {
                a.id = core.clone();
                b.id = core;
            }
        }
        io::write_records(&args.io.output, SeqFormat::Fastq, &args.io.compression, &r1)?;
        io::write_records(out2, SeqFormat::Fastq, &args.io.compression, &r2)?;
        return Ok(());
    }

    let in_path = args.io.input.as_deref();
    let fmt = SeqFormat::from_arg(&args.io.format).unwrap_or(SeqFormat::detect(in_path)?);
    let mut recs = io::read_records(in_path, fmt, &args.io.compression)?;
    for (i, r) in recs.iter_mut().enumerate() {
        r.id = rename_id(&args, mode, i, &r.id, &mapping, regex.as_ref(), replacement)?;
    }
    io::write_records(&args.io.output, fmt, &args.io.compression, &recs)
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
    let txt = fs::read_to_string(path)
        .with_context(|| format!("failed to read mapping file '{path}'"))?;
    let mut map = HashMap::new();
    for (i, raw_line) in txt.lines().enumerate() {
        let line_no = i + 1;
        let line = raw_line.trim_end();
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
    mapping: &Option<HashMap<String, String>>,
    regex: Option<&Regex>,
    replacement: Option<&str>,
) -> Result<String> {
    match mode {
        RenameMode::Generate => Ok(format_name(args, args.start + idx)),
        RenameMode::Map => Ok(mapping
            .as_ref()
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
