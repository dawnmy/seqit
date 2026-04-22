use anyhow::Result;

use crate::{
    cli::{SortArgs, SortBy},
    formats::{SeqFormat, SeqRecord},
    io,
};

pub fn run(args: SortArgs) -> Result<()> {
    let in_path = args.io.input.as_deref();
    let fmt = SeqFormat::from_arg(&args.io.format).unwrap_or(SeqFormat::detect(in_path)?);
    let mut recs = io::read_records(in_path, fmt, &args.io.compression)?;
    recs.sort_by(|a, b| key(a, &args.by).cmp(&key(b, &args.by)));
    if args.reverse {
        recs.reverse();
    }
    io::write_records(&args.io.output, fmt, &args.io.compression, &recs)
}

fn key<'a>(r: &'a SeqRecord, by: &SortBy) -> Vec<u8> {
    match by {
        SortBy::Id | SortBy::Name => r.id.as_bytes().to_vec(),
        SortBy::Len => format!("{:020}", r.seq.len()).into_bytes(),
        SortBy::Seq => r.seq.clone(),
    }
}
