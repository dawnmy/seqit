use std::collections::HashSet;
use std::io::Write;

use anyhow::{bail, Result};

use crate::cli::SeqArgs;
use crate::formats::{SeqFormat, SeqRecord};
use crate::io;

pub fn run(args: SeqArgs) -> Result<()> {
    let in_path = args.io.input.as_deref();
    let (fmt, br) = io::open_seq_reader(in_path, &args.io.format, &args.io.compression)?;
    if !matches!(fmt, SeqFormat::Fasta | SeqFormat::Fastq) {
        bail!("seq currently supports FASTA/FASTQ input");
    }
    if (args.min_qual >= 0.0 || args.max_qual >= 0.0) && !matches!(fmt, SeqFormat::Fastq) {
        bail!("--min-qual/--max-qual require FASTQ records with qualities");
    }

    let w = io::open_writer(&args.io.output)?;
    let w = io::wrap_compress(w, &args.io.output, &args.io.compression)?;
    let mut bw = io::buffered_writer(w);

    let gap_letters = args
        .remove_gaps
        .then(|| parse_gap_letters(&args.gap_letters));
    let needs_owned = needs_owned_record(&args);
    let parse_headers = needs_parsed_header(&args);
    let plain_output = can_write_plain_records(&args);
    io::for_each_record_from_reader_with_header_parsing(br, fmt, parse_headers, |rec| {
        if plain_output {
            io::write_record_ref(&mut bw, fmt, &rec)?;
            return Ok(());
        }

        if needs_owned {
            let mut rec = rec.to_owned_record()?;
            if process_owned_record(&args, &gap_letters, &mut rec)? {
                write_record(&args, fmt, &mut bw, &rec)?;
            }
            return Ok(());
        }

        if process_record_ref(&args, &rec)? {
            write_record_ref(&args, fmt, &mut bw, &rec)?;
        }
        Ok(())
    })?;

    bw.flush()?;
    Ok(())
}

fn needs_owned_record(args: &SeqArgs) -> bool {
    args.rev || args.comp || args.revcomp || args.upper || args.lower || args.remove_gaps
}

fn needs_parsed_header(args: &SeqArgs) -> bool {
    args.only_id || (args.name && !args.full_name) || args.validate_seq
}

fn can_write_plain_records(args: &SeqArgs) -> bool {
    !needs_owned_record(args)
        && !needs_parsed_header(args)
        && args.min_len.is_none()
        && args.max_len.is_none()
        && args.min_qual < 0.0
        && args.max_qual < 0.0
        && !args.seq
        && !args.name
        && !args.color
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

fn process_record_ref(args: &SeqArgs, rec: &io::SeqRecordRef<'_>) -> Result<bool> {
    let len = rec.seq.len();
    let len_ok = args.min_len.map(|m| len >= m).unwrap_or(true)
        && args.max_len.map(|m| len <= m).unwrap_or(true);
    if !len_ok {
        return Ok(false);
    }

    let qual_filter_enabled = args.min_qual >= 0.0 || args.max_qual >= 0.0;
    if qual_filter_enabled {
        let Some(q) = rec.qual else {
            bail!("--min-qual/--max-qual require FASTQ records with qualities");
        };
        let avg = average_quality(q, args.qual_ascii_base);
        if (args.min_qual >= 0.0 && avg < args.min_qual)
            || (args.max_qual >= 0.0 && avg > args.max_qual)
        {
            return Ok(false);
        }
    }

    if args.validate_seq {
        if let Some(ch) = find_invalid_iupac(&rec.seq) {
            bail!(
                "record '{}' contains invalid base '{}'",
                rec.id_str()?,
                ch as char
            );
        }
    }

    Ok(true)
}

fn process_owned_record(
    args: &SeqArgs,
    gap_letters: &Option<HashSet<u8>>,
    rec: &mut SeqRecord,
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
        crate::utils::reverse_complement_in_place(&mut rec.seq);
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
            crate::utils::complement_in_place(&mut rec.seq);
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
    rec: &SeqRecord,
) -> Result<()> {
    if args.seq {
        write_seq_line(bw, &rec.seq, args.color)?;
        return Ok(());
    }
    if args.only_id {
        write_line_bytes(bw, rec.id.as_bytes())?;
        return Ok(());
    }
    if args.name {
        if args.full_name {
            write_full_name_line(bw, rec)?;
        } else {
            write_line_bytes(bw, rec.id.as_bytes())?;
        }
        return Ok(());
    }

    match fmt {
        SeqFormat::Fasta => {
            bw.write_all(b">")?;
            write_full_name_line(bw, rec)?;
            write_seq_line(bw, &rec.seq, args.color)?;
        }
        SeqFormat::Fastq => {
            bw.write_all(b"@")?;
            write_full_name_line(bw, rec)?;
            write_seq_line(bw, &rec.seq, args.color)?;
            bw.write_all(b"+\n")?;
            let Some(q) = &rec.qual else {
                bail!("record '{}' missing qualities for FASTQ output", rec.id);
            };
            write_line_bytes(bw, q)?;
        }
        _ => bail!("format {:?} not yet implemented for seq output", fmt),
    }
    Ok(())
}

fn write_record_ref(
    args: &SeqArgs,
    fmt: SeqFormat,
    bw: &mut impl Write,
    rec: &io::SeqRecordRef<'_>,
) -> Result<()> {
    if args.seq {
        write_seq_line(bw, &rec.seq, args.color)?;
        return Ok(());
    }
    if args.only_id {
        write_line_bytes(bw, rec.id)?;
        return Ok(());
    }
    if args.name {
        if args.full_name {
            io::write_header_bytes(bw, rec.id, rec.desc)?;
            bw.write_all(b"\n")?;
        } else {
            write_line_bytes(bw, rec.id)?;
        }
        return Ok(());
    }

    io::write_record_ref(bw, fmt, rec)
}

fn write_line_bytes(writer: &mut impl Write, line: &[u8]) -> Result<()> {
    writer.write_all(line)?;
    writer.write_all(b"\n")?;
    Ok(())
}

fn write_full_name_line(writer: &mut impl Write, rec: &SeqRecord) -> Result<()> {
    writer.write_all(rec.id.as_bytes())?;
    if let Some(desc) = &rec.desc {
        writer.write_all(b" ")?;
        writer.write_all(desc.as_bytes())?;
    }
    writer.write_all(b"\n")?;
    Ok(())
}

fn write_seq_line(writer: &mut impl Write, seq: &[u8], color: bool) -> Result<()> {
    if !color {
        return write_line_bytes(writer, seq);
    }
    write_line_bytes(writer, render_colored_seq(seq).as_bytes())
}

fn render_colored_seq(seq: &[u8]) -> String {
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
