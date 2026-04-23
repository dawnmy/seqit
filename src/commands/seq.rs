use std::collections::HashSet;
use std::io::{BufReader, BufWriter, Write};

use anyhow::{bail, Result};
use bio::io::{fasta, fastq};

use crate::cli::{FormatArg, SeqArgs};
use crate::formats::SeqFormat;
use crate::io;

pub fn run(args: SeqArgs) -> Result<()> {
    let in_path = args.io.input.as_deref();
    if matches!(args.io.format, FormatArg::Auto) && matches!(in_path, None | Some("-")) {
        let (fmt, recs) =
            io::read_records_with_format(in_path, &args.io.format, &args.io.compression)?;
        return run_from_records(args, fmt, recs);
    }

    let fmt = io::resolve_seq_format(in_path, &args.io.format, &args.io.compression)?;
    if !matches!(fmt, SeqFormat::Fasta | SeqFormat::Fastq) {
        bail!("seq currently supports FASTA/FASTQ input");
    }
    if (args.min_qual >= 0.0 || args.max_qual >= 0.0) && !matches!(fmt, SeqFormat::Fastq) {
        bail!("--min-qual/--max-qual require FASTQ records with qualities");
    }

    let r = io::open_reader(in_path)?;
    let r = io::wrap_decompress(r, in_path, &args.io.compression)?;
    let br = BufReader::new(r);
    let w = io::open_writer(&args.io.output)?;
    let w = io::wrap_compress(w, &args.io.output, &args.io.compression)?;
    let mut bw = BufWriter::new(w);

    let gap_letters = args
        .remove_gaps
        .then(|| parse_gap_letters(&args.gap_letters));
    match fmt {
        SeqFormat::Fasta => {
            let fa = fasta::Reader::new(br);
            for r in fa.records() {
                let r = r?;
                let mut rec = crate::formats::SeqRecord {
                    id: r.id().to_string(),
                    desc: r.desc().map(|d| d.to_string()),
                    seq: r.seq().to_vec(),
                    qual: None,
                };
                if !process_record(&args, &gap_letters, &mut rec)? {
                    continue;
                }
                write_record(&args, fmt, &mut bw, &rec)?;
            }
        }
        SeqFormat::Fastq => {
            let fq = fastq::Reader::new(br);
            for r in fq.records() {
                let r = r?;
                let mut rec = crate::formats::SeqRecord {
                    id: r.id().to_string(),
                    desc: r.desc().map(|d| d.to_string()),
                    seq: r.seq().to_vec(),
                    qual: Some(r.qual().to_vec()),
                };
                if !process_record(&args, &gap_letters, &mut rec)? {
                    continue;
                }
                write_record(&args, fmt, &mut bw, &rec)?;
            }
        }
        _ => bail!("format {:?} not yet implemented for seq output", fmt),
    }

    bw.flush()?;
    Ok(())
}

fn run_from_records(
    args: SeqArgs,
    fmt: SeqFormat,
    recs: Vec<crate::formats::SeqRecord>,
) -> Result<()> {
    if !matches!(fmt, SeqFormat::Fasta | SeqFormat::Fastq) {
        bail!("seq currently supports FASTA/FASTQ input");
    }
    if (args.min_qual >= 0.0 || args.max_qual >= 0.0) && !matches!(fmt, SeqFormat::Fastq) {
        bail!("--min-qual/--max-qual require FASTQ records with qualities");
    }
    let w = io::open_writer(&args.io.output)?;
    let w = io::wrap_compress(w, &args.io.output, &args.io.compression)?;
    let mut bw = BufWriter::new(w);
    let gap_letters = args
        .remove_gaps
        .then(|| parse_gap_letters(&args.gap_letters));
    for mut rec in recs {
        if !process_record(&args, &gap_letters, &mut rec)? {
            continue;
        }
        write_record(&args, fmt, &mut bw, &rec)?;
    }
    bw.flush()?;
    Ok(())
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

fn process_record(
    args: &SeqArgs,
    gap_letters: &Option<HashSet<u8>>,
    rec: &mut crate::formats::SeqRecord,
) -> Result<bool> {
    let len = rec.seq.len();
    let len_ok = args.min_len.map(|m| len >= m).unwrap_or(true)
        && args.max_len.map(|m| len <= m).unwrap_or(true);
    if !len_ok {
        return Ok(false);
    }

    let qual_filter_enabled = args.min_qual >= 0.0 || args.max_qual >= 0.0;
    if qual_filter_enabled {
        let Some(q) = &rec.qual else {
            bail!("--min-qual/--max-qual require FASTQ records with qualities");
        };
        let avg = average_quality(q, args.qual_ascii_base);
        if (args.min_qual >= 0.0 && avg < args.min_qual)
            || (args.max_qual >= 0.0 && avg > args.max_qual)
        {
            return Ok(false);
        }
    }

    if args.revcomp {
        rec.seq = crate::utils::revcomp(&rec.seq);
        if let Some(q) = &mut rec.qual {
            q.reverse();
        }
    } else {
        if args.rev {
            rec.seq.reverse();
            if let Some(q) = &mut rec.qual {
                q.reverse();
            }
        }
        if args.comp {
            rec.seq = crate::utils::revcomp(&rec.seq).into_iter().rev().collect();
        }
    }
    if args.upper {
        rec.seq.make_ascii_uppercase();
    }
    if args.lower {
        rec.seq.make_ascii_lowercase();
    }
    if let Some(gaps) = gap_letters {
        rec.seq.retain(|b| !gaps.contains(b));
    }
    if args.validate_seq {
        if let Some(ch) = find_invalid_iupac(&rec.seq) {
            bail!("record '{}' contains invalid base '{}'", rec.id, ch as char);
        }
    }

    Ok(true)
}

fn write_record(
    args: &SeqArgs,
    fmt: SeqFormat,
    bw: &mut impl Write,
    rec: &crate::formats::SeqRecord,
) -> Result<()> {
    if args.seq {
        write_line(bw, &render_seq(&rec.seq, args.color))?;
        return Ok(());
    }
    if args.only_id {
        write_line(bw, &rec.id)?;
        return Ok(());
    }
    if args.name {
        if args.full_name {
            write_line(bw, &rec.name())?;
        } else {
            write_line(bw, &rec.id)?;
        }
        return Ok(());
    }

    match fmt {
        SeqFormat::Fasta => {
            write_line(bw, &format!(">{}", rec.name()))?;
            write_line(bw, &render_seq(&rec.seq, args.color))?;
        }
        SeqFormat::Fastq => {
            write_line(bw, &format!("@{}", rec.name()))?;
            write_line(bw, &render_seq(&rec.seq, args.color))?;
            write_line(bw, "+")?;
            let Some(q) = &rec.qual else {
                bail!("record '{}' missing qualities for FASTQ output", rec.id);
            };
            write_line(bw, &String::from_utf8_lossy(q))?;
        }
        _ => bail!("format {:?} not yet implemented for seq output", fmt),
    }
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
