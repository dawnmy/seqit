use anyhow::{bail, Context, Result};

use crate::{cli::HeadArgs, formats::SeqFormat, io, utils};

pub fn run(args: HeadArgs) -> Result<()> {
    let n = resolve_take(args.num, args.proportion)?;
    utils::validate_input_mode(
        "head",
        args.io.input.as_deref(),
        args.input1.as_deref(),
        args.input2.as_deref(),
    )?;
    if let (Some(in1), Some(in2)) = (args.input1.as_deref(), args.input2.as_deref()) {
        let out2 = args
            .output2
            .as_deref()
            .context("paired head requires --output2")?;
        if let Some(num) = args.num {
            return run_paired_streaming_num(&args, in1, in2, out2, num);
        }
        if in1 == "-" || in2 == "-" {
            bail!("head --proportion for paired stdin is not supported; use --num for streaming paired input");
        }
        let total = count_valid_pairs(&args, in1, in2)?;
        return run_paired_streaming_num(&args, in1, in2, out2, n(total));
    }

    let in_path = args.io.input.as_deref();
    if let Some(num) = args.num {
        return run_single_streaming_num(&args, in_path, num);
    }
    if !matches!(in_path, None | Some("-")) {
        let (fmt, r) = io::open_seq_reader(in_path, &args.io.format, &args.io.compression)?;
        let total = count_records(fmt, r)?;
        return run_single_streaming_num(&args, in_path, n(total));
    }
    let (fmt, recs) = io::read_records_with_format(in_path, &args.io.format, &args.io.compression)?;
    let keep = n(recs.len());
    io::write_records(&args.io.output, fmt, &args.io.compression, &recs[..keep])
}

fn run_single_streaming_num(args: &HeadArgs, in_path: Option<&str>, n: usize) -> Result<()> {
    let (fmt, r) = io::open_seq_reader(in_path, &args.io.format, &args.io.compression)?;
    let w = io::open_writer(&args.io.output)?;
    let w = io::wrap_compress(w, &args.io.output, &args.io.compression)?;
    let mut out = io::buffered_writer(w);
    let mut kept = 0usize;
    io::for_each_record_from_reader_until(r, fmt, |rec| {
        if kept >= n {
            return Ok(false);
        }
        io::write_record_ref(&mut out, fmt, &rec)?;
        kept += 1;
        Ok(kept < n)
    })?;
    Ok(())
}

fn run_paired_streaming_num(
    args: &HeadArgs,
    in1: &str,
    in2: &str,
    out2: &str,
    n: usize,
) -> Result<()> {
    let w1 = io::open_writer(&args.io.output)?;
    let w1 = io::wrap_compress_for_streams(w1, &args.io.output, &args.io.compression, 2)?;
    let w2 = io::open_writer(out2)?;
    let w2 = io::wrap_compress_for_streams(w2, out2, &args.io.compression, 2)?;
    let mut w1 = io::buffered_writer(w1);
    let mut w2 = io::buffered_writer(w2);

    let mut kept = 0usize;
    let mut invalid = 0usize;
    if n > 0 {
        io::for_each_fastq_pair(in1, in2, &args.io.compression, |pair| {
            if kept >= n {
                return Ok(false);
            }
            match pair {
                io::FastqPair::Both(a, b) => {
                    let a_id = a.id_str()?;
                    let b_id = b.id_str()?;
                    if pair_key(a_id) != pair_key(b_id) {
                        invalid += 1;
                        if !args.allow_unpaired {
                            bail!(
                                "paired inputs contain invalid/unpaired records (first mismatch: R1={} R2={}); rerun with --allow-unpaired to continue",
                                a_id,
                                b_id
                            );
                        }
                        return Ok(true);
                    }
                    io::write_record_ref(&mut w1, SeqFormat::Fastq, &a)?;
                    io::write_record_ref(&mut w2, SeqFormat::Fastq, &b)?;
                    kept += 1;
                    Ok(kept < n)
                }
                io::FastqPair::Left(a) => {
                    invalid += 1;
                    if !args.allow_unpaired {
                        bail!(
                            "paired inputs contain unmatched record in R1 ({}); rerun with --allow-unpaired to continue",
                            a.id_str()?
                        );
                    }
                    Ok(true)
                }
                io::FastqPair::Right(b) => {
                    invalid += 1;
                    if !args.allow_unpaired {
                        bail!(
                            "paired inputs contain unmatched record in R2 ({}); rerun with --allow-unpaired to continue",
                            b.id_str()?
                        );
                    }
                    Ok(true)
                }
            }
        })?;
    }
    if invalid > 0 && args.allow_unpaired {
        eprintln!("warning: skipped {invalid} invalid/unpaired paired records");
    }
    Ok(())
}

fn pair_key(id: &str) -> &str {
    id.strip_suffix("/1")
        .or_else(|| id.strip_suffix("/2"))
        .unwrap_or(id)
}

fn count_records(fmt: SeqFormat, br: io::DynBufReader) -> Result<usize> {
    let mut total = 0usize;
    io::for_each_record_from_reader(br, fmt, |_| {
        total += 1;
        Ok(())
    })?;
    Ok(total)
}

fn count_valid_pairs(args: &HeadArgs, in1: &str, in2: &str) -> Result<usize> {
    let mut valid = 0usize;
    io::for_each_fastq_pair(in1, in2, &args.io.compression, |pair| {
        match pair {
            io::FastqPair::Both(a, b) => {
                let a_id = a.id_str()?;
                let b_id = b.id_str()?;
                if pair_key(a_id) == pair_key(b_id) {
                    valid += 1;
                } else if !args.allow_unpaired {
                    bail!(
                        "paired inputs contain invalid/unpaired records (first mismatch: R1={} R2={}); rerun with --allow-unpaired to continue",
                        a_id,
                        b_id
                    );
                }
            }
            io::FastqPair::Left(a) if !args.allow_unpaired => {
                bail!(
                    "paired inputs contain unmatched record in R1 ({}); rerun with --allow-unpaired to continue",
                    a.id_str()?
                );
            }
            io::FastqPair::Right(b) if !args.allow_unpaired => {
                bail!(
                    "paired inputs contain unmatched record in R2 ({}); rerun with --allow-unpaired to continue",
                    b.id_str()?
                );
            }
            _ => {}
        }
        Ok(true)
    })?;
    Ok(valid)
}

fn resolve_take(
    num: Option<usize>,
    proportion: Option<f64>,
) -> Result<Box<dyn Fn(usize) -> usize>> {
    match (num, proportion) {
        (Some(_), Some(_)) => bail!("choose only one of --num and --proportion"),
        (None, None) => bail!("one of --num or --proportion is required"),
        (Some(n), None) => Ok(Box::new(move |len| n.min(len))),
        (None, Some(p)) => {
            if !(0.0..=1.0).contains(&p) {
                bail!("--proportion must be in [0,1]");
            }
            Ok(Box::new(move |len| ((len as f64) * p).ceil() as usize))
        }
    }
}
