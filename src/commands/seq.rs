use std::collections::HashSet;
use std::io::{BufWriter, Write};

use anyhow::{bail, Result};
use rayon::prelude::*;

use crate::cli::SeqArgs;
use crate::formats::SeqFormat;
use crate::io;

pub fn run(args: SeqArgs) -> Result<()> {
    let in_path = args.io.input.as_deref();
    let fmt = SeqFormat::from_arg(&args.io.format).unwrap_or(SeqFormat::detect(in_path)?);
    if !matches!(fmt, SeqFormat::Fasta | SeqFormat::Fastq) {
        bail!("seq currently supports FASTA/FASTQ input");
    }
    let mut recs = io::read_records(in_path, fmt, &args.io.compression)?;
    let qual_filter_enabled = args.min_qual >= 0.0 || args.max_qual >= 0.0;
    if qual_filter_enabled && recs.iter().any(|r| r.qual.is_none()) {
        bail!("--min-qual/--max-qual require FASTQ records with qualities");
    }

    recs = recs
        .into_par_iter()
        .filter(|r| {
            let len = r.seq.len();
            let len_ok = args.min_len.map(|m| len >= m).unwrap_or(true)
                && args.max_len.map(|m| len <= m).unwrap_or(true);
            let qual_ok = if qual_filter_enabled {
                r.qual
                    .as_ref()
                    .map(|q| {
                        let avg = average_quality(q, args.qual_ascii_base);
                        (args.min_qual < 0.0 || avg >= args.min_qual)
                            && (args.max_qual < 0.0 || avg <= args.max_qual)
                    })
                    .unwrap_or(false)
            } else {
                true
            };
            len_ok && qual_ok
        })
        .collect();

    let gap_letters = if args.remove_gaps {
        Some(parse_gap_letters(&args.gap_letters))
    } else {
        None
    };

    recs.par_iter_mut().for_each(|r| {
        if args.revcomp {
            r.seq = crate::utils::revcomp(&r.seq);
            if let Some(q) = &mut r.qual {
                q.reverse();
            }
        } else {
            if args.rev {
                r.seq.reverse();
                if let Some(q) = &mut r.qual {
                    q.reverse();
                }
            }
            if args.comp {
                r.seq = crate::utils::revcomp(&r.seq).into_iter().rev().collect();
            }
        }
        if args.upper {
            r.seq.make_ascii_uppercase();
        }
        if args.lower {
            r.seq.make_ascii_lowercase();
        }
        if let Some(gaps) = &gap_letters {
            r.seq.retain(|b| !gaps.contains(b));
        }
    });

    if args.validate_seq {
        for r in &recs {
            if let Some(ch) = find_invalid_iupac(&r.seq) {
                bail!("record '{}' contains invalid base '{}'", r.id, ch as char);
            }
        }
    }

    write_output(&args, fmt, &recs)
}

fn average_quality(qual: &[u8], ascii_base: u8) -> f64 {
    if qual.is_empty() {
        return 0.0;
    }
    let sum: i64 = qual
        .iter()
        .map(|q| i64::from(*q) - i64::from(ascii_base))
        .sum();
    sum as f64 / qual.len() as f64
}

fn parse_gap_letters(gap_letters: &str) -> HashSet<u8> {
    gap_letters.as_bytes().iter().copied().collect()
}

fn find_invalid_iupac(seq: &[u8]) -> Option<u8> {
    seq.iter().copied().find(|b| {
        !matches!(
            b.to_ascii_uppercase(),
            b'A' | b'C'
                | b'G'
                | b'T'
                | b'U'
                | b'R'
                | b'Y'
                | b'S'
                | b'W'
                | b'K'
                | b'M'
                | b'B'
                | b'D'
                | b'H'
                | b'V'
                | b'N'
        )
    })
}

fn write_output(args: &SeqArgs, fmt: SeqFormat, recs: &[crate::formats::SeqRecord]) -> Result<()> {
    let w = io::open_writer(&args.io.output)?;
    let w = io::wrap_compress(w, &args.io.output, &args.io.compression)?;
    let mut bw = BufWriter::new(w);

    for r in recs {
        if args.seq {
            write_line(&mut bw, &render_seq(&r.seq, args.color))?;
            continue;
        }
        if args.only_id {
            write_line(&mut bw, &r.id)?;
            continue;
        }
        if args.name {
            if args.full_name {
                write_line(&mut bw, &r.name())?;
            } else {
                write_line(&mut bw, &r.id)?;
            }
            continue;
        }

        match fmt {
            SeqFormat::Fasta => {
                write_line(&mut bw, &format!(">{}", r.name()))?;
                write_line(&mut bw, &render_seq(&r.seq, args.color))?;
            }
            SeqFormat::Fastq => {
                write_line(&mut bw, &format!("@{}", r.name()))?;
                write_line(&mut bw, &render_seq(&r.seq, args.color))?;
                write_line(&mut bw, "+")?;
                let Some(q) = &r.qual else {
                    bail!("record '{}' missing qualities for FASTQ output", r.id);
                };
                write_line(&mut bw, &String::from_utf8_lossy(q))?;
            }
            _ => bail!("format {:?} not yet implemented for seq output", fmt),
        }
    }

    bw.flush()?;
    Ok(())
}

fn write_line(writer: &mut impl Write, line: &str) -> Result<()> {
    writer.write_all(line.as_bytes())?;
    writer.write_all(b"\n")?;
    Ok(())
}

fn render_seq(seq: &[u8], color: bool) -> String {
    if !color {
        return String::from_utf8_lossy(seq).into_owned();
    }
    let mut out = String::with_capacity(seq.len() * 8);
    for b in seq {
        let c = *b as char;
        let code = match b.to_ascii_uppercase() {
            b'A' => "32",
            b'C' => "34",
            b'G' => "33",
            b'T' | b'U' => "31",
            _ => "37",
        };
        out.push_str("\u{1b}[");
        out.push_str(code);
        out.push('m');
        out.push(c);
        out.push_str("\u{1b}[0m");
    }
    out
}
