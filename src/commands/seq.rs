use anyhow::Result;
use rayon::prelude::*;

use crate::cli::SeqArgs;
use crate::formats::SeqFormat;
use crate::io;

pub fn run(args: SeqArgs) -> Result<()> {
    let in_path = args.io.input.as_deref();
    let fmt = SeqFormat::from_arg(&args.io.format).unwrap_or(SeqFormat::detect(in_path)?);
    let mut recs = io::read_records(in_path, fmt, &args.io.compression)?;

    recs = recs
        .into_par_iter()
        .filter(|r| {
            let len = r.seq.len();
            args.min_len.map(|m| len >= m).unwrap_or(true)
                && args.max_len.map(|m| len <= m).unwrap_or(true)
        })
        .collect();

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
    });

    io::write_records(&args.io.output, fmt, &args.io.compression, &recs)
}
