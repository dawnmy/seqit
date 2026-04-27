use anyhow::{bail, Context, Result};
use std::collections::VecDeque;

use crate::{
    cli::TailArgs,
    formats::{SeqFormat, SeqRecord},
    io, utils,
};

type PairParts = (SeqRecord, SeqRecord);

pub fn run(args: TailArgs) -> Result<()> {
    let _ = resolve_take(args.num, args.proportion)?;
    utils::validate_input_mode(
        "tail",
        args.io.input.as_deref(),
        args.input1.as_deref(),
        args.input2.as_deref(),
    )?;
    if let (Some(in1), Some(in2)) = (args.input1.as_deref(), args.input2.as_deref()) {
        let out2 = args
            .output2
            .as_deref()
            .context("paired tail requires --output2")?;
        return run_paired_tail_streaming(&args, in1, in2, out2);
    }

    let in_path = args.io.input.as_deref();
    run_single_tail_streaming(&args, in_path)
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

fn run_single_tail_streaming(args: &TailArgs, in_path: Option<&str>) -> Result<()> {
    let mode = TailMode::from_args(args)?;
    let (fmt, br) = io::open_seq_reader(in_path, &args.io.format, &args.io.compression)?;
    let w = io::open_writer(&args.io.output)?;
    let w = io::wrap_compress(w, &args.io.output, &args.io.compression)?;
    let mut bw = io::buffered_writer(w);

    if !matches!(fmt, SeqFormat::Fasta | SeqFormat::Fastq) {
        bail!("tail currently supports FASTA/FASTQ input");
    }
    let mut ring: VecDeque<SeqRecord> = VecDeque::with_capacity(mode.initial_capacity());
    let mut seen = 0usize;
    io::for_each_record_from_reader(br, fmt, |rec| {
        seen += 1;
        push_tail(&mut ring, mode, seen, rec.to_owned_record()?);
        Ok(())
    })?;
    for rec in ring {
        io::write_record(&mut bw, fmt, &rec)?;
    }
    Ok(())
}

fn run_paired_tail_streaming(args: &TailArgs, in1: &str, in2: &str, out2: &str) -> Result<()> {
    let mode = TailMode::from_args(args)?;
    let mut ring: VecDeque<PairParts> = VecDeque::with_capacity(mode.initial_capacity());
    let mut valid_seen = 0usize;
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
                valid_seen += 1;
                push_tail(
                    &mut ring,
                    mode,
                    valid_seen,
                    (a.to_owned_record()?, b.to_owned_record()?),
                );
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
    report_invalid_pairs(invalid_count, invalid_preview, args.allow_unpaired)?;

    let w1 = io::open_writer(&args.io.output)?;
    let w1 = io::wrap_compress_for_streams(w1, &args.io.output, &args.io.compression, 2)?;
    let w2 = io::open_writer(out2)?;
    let w2 = io::wrap_compress_for_streams(w2, out2, &args.io.compression, 2)?;
    let mut w1 = io::buffered_writer(w1);
    let mut w2 = io::buffered_writer(w2);
    for (a, b) in ring {
        io::write_record(&mut w1, SeqFormat::Fastq, &a)?;
        io::write_record(&mut w2, SeqFormat::Fastq, &b)?;
    }
    Ok(())
}

fn pair_key(id: &str) -> &str {
    id.strip_suffix("/1")
        .or_else(|| id.strip_suffix("/2"))
        .unwrap_or(id)
}

#[derive(Clone, Copy)]
enum TailMode {
    Fixed(usize),
    Proportion(f64),
}

impl TailMode {
    fn from_args(args: &TailArgs) -> Result<Self> {
        match (args.num, args.proportion) {
            (Some(n), None) => Ok(Self::Fixed(n)),
            (None, Some(p)) => Ok(Self::Proportion(p)),
            _ => unreachable!("tail arguments are validated before streaming"),
        }
    }

    fn initial_capacity(self) -> usize {
        match self {
            Self::Fixed(n) => n,
            Self::Proportion(_) => 0,
        }
    }

    fn keep_after_seen(self, seen: usize) -> usize {
        match self {
            Self::Fixed(n) => n,
            Self::Proportion(p) => ((seen as f64) * p).ceil() as usize,
        }
    }
}

fn push_tail<T>(ring: &mut VecDeque<T>, mode: TailMode, seen: usize, item: T) {
    let keep = mode.keep_after_seen(seen);
    if keep == 0 {
        return;
    }
    while ring.len() >= keep {
        ring.pop_front();
    }
    ring.push_back(item);
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
