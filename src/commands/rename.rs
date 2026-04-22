use anyhow::{Context, Result};

use crate::{cli::RenameArgs, formats::SeqFormat, io, pairs};

pub fn run(args: RenameArgs) -> Result<()> {
    if let (Some(in1), Some(in2)) = (args.in1.as_deref(), args.in2.as_deref()) {
        let mut r1 = io::read_records(Some(in1), SeqFormat::Fastq, &args.io.compression)?;
        let mut r2 = io::read_records(Some(in2), SeqFormat::Fastq, &args.io.compression)?;
        pairs::validate_pair_counts(&r1, &r2)?;
        let out2 = args
            .output2
            .as_deref()
            .context("paired rename requires --output2")?;
        for (i, (a, b)) in r1.iter_mut().zip(r2.iter_mut()).enumerate() {
            let n = args.start + i;
            let core = format_name(&args, n);
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
        let n = args.start + i;
        r.id = format_name(&args, n);
    }
    io::write_records(&args.io.output, fmt, &args.io.compression, &recs)
}

fn format_name(args: &RenameArgs, n: usize) -> String {
    if let Some(tpl) = &args.template {
        tpl.replace("{prefix}", &args.prefix)
            .replace("{n}", &format!("{:0width$}", n, width = args.width))
    } else {
        format!("{}{:0width$}", args.prefix, n, width = args.width)
    }
}
