use anyhow::{bail, Context, Result};

use crate::{cli::HeadArgs, formats::SeqFormat, io, pairs};

pub fn run(args: HeadArgs) -> Result<()> {
    let n = resolve_take(args.num, args.proportion)?;
    if let (Some(in1), Some(in2)) = (args.input1.as_deref(), args.input2.as_deref()) {
        let r1 = io::read_records(Some(in1), SeqFormat::Fastq, &args.io.compression)?;
        let r2 = io::read_records(Some(in2), SeqFormat::Fastq, &args.io.compression)?;
        let (r1, r2) = pairs::prepare_paired_records(r1, r2, args.allow_unpaired)?;
        let keep = n(r1.len());
        let out2 = args
            .output2
            .as_deref()
            .context("paired head requires --output2")?;
        io::write_records(
            &args.io.output,
            SeqFormat::Fastq,
            &args.io.compression,
            &r1[..keep],
        )?;
        io::write_records(out2, SeqFormat::Fastq, &args.io.compression, &r2[..keep])?;
        return Ok(());
    }

    let in_path = args.io.input.as_deref();
    let fmt = SeqFormat::from_arg(&args.io.format).unwrap_or(SeqFormat::detect(in_path)?);
    let recs = io::read_records(in_path, fmt, &args.io.compression)?;
    let keep = n(recs.len());
    io::write_records(&args.io.output, fmt, &args.io.compression, &recs[..keep])
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
