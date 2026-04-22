use std::cmp::Ordering;

use anyhow::Result;

use crate::{
    cli::{SortArgs, SortBy},
    formats::{SeqFormat, SeqRecord},
    io,
    utils::parse_mem_bytes,
};

pub fn run(args: SortArgs) -> Result<()> {
    let in_path = args.io.input.as_deref();
    let fmt = SeqFormat::from_arg(&args.io.format).unwrap_or(SeqFormat::detect(in_path)?);
    let _mem_limit_bytes = parse_mem_bytes(&args.mem);
    let mut recs = io::read_records(in_path, fmt, &args.io.compression)?;
    recs.sort_by(|a, b| compare_records(a, b, &args.by, args.numeric));
    if args.reverse {
        recs.reverse();
    }
    io::write_records(&args.io.output, fmt, &args.io.compression, &recs)
}

fn key(r: &SeqRecord, by: &SortBy) -> Vec<u8> {
    match by {
        SortBy::Id => r.id.as_bytes().to_vec(),
        SortBy::Name => r.name().into_bytes(),
        SortBy::Len => format!("{:020}", r.seq.len()).into_bytes(),
        SortBy::Seq => r.seq.clone(),
    }
}

fn compare_records(a: &SeqRecord, b: &SeqRecord, by: &SortBy, numeric: bool) -> Ordering {
    if !numeric || matches!(by, SortBy::Seq) {
        return key(a, by).cmp(&key(b, by));
    }
    match by {
        SortBy::Len => a.seq.len().cmp(&b.seq.len()),
        SortBy::Id | SortBy::Name => {
            let ak = key(a, by);
            let bk = key(b, by);
            match (
                String::from_utf8_lossy(&ak).parse::<u128>(),
                String::from_utf8_lossy(&bk).parse::<u128>(),
            ) {
                (Ok(x), Ok(y)) => x.cmp(&y),
                _ => ak.cmp(&bk),
            }
        }
        SortBy::Seq => Ordering::Equal,
    }
}
