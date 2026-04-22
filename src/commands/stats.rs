use anyhow::{bail, Result};
use serde::Serialize;

use crate::cli::StatsArgs;
use crate::formats::SeqFormat;
use crate::io;

#[derive(Debug, Serialize, Clone)]
struct StatRow {
    file: String,
    format: String,
    records: usize,
    total_bases: usize,
    min_len: usize,
    max_len: usize,
    mean_len: f64,
    gc_pct: f64,
    n50: usize,
}

pub fn run(args: StatsArgs) -> Result<()> {
    let inputs = if args.inputs.is_empty() {
        vec!["-".to_string()]
    } else {
        args.inputs.clone()
    };
    let mut rows = Vec::new();
    for p in inputs {
        let fmt = SeqFormat::from_arg(&args.format).unwrap_or(SeqFormat::detect(Some(&p))?);
        match fmt {
            SeqFormat::Fasta | SeqFormat::Fastq => {
                let recs = io::read_records(Some(&p), fmt, &crate::cli::CompressionArg::Auto)?;
                let lens: Vec<usize> = recs.iter().map(|r| r.seq.len()).collect();
                let records = lens.len();
                let total_bases: usize = lens.iter().sum();
                let min_len = *lens.iter().min().unwrap_or(&0);
                let max_len = *lens.iter().max().unwrap_or(&0);
                let mean_len = if records > 0 {
                    total_bases as f64 / records as f64
                } else {
                    0.0
                };
                let gc = recs
                    .iter()
                    .map(|r| {
                        bytecount::count(&r.seq, b'G')
                            + bytecount::count(&r.seq, b'g')
                            + bytecount::count(&r.seq, b'C')
                            + bytecount::count(&r.seq, b'c')
                    })
                    .sum::<usize>();
                let gc_pct = if total_bases > 0 {
                    (gc as f64) * 100.0 / (total_bases as f64)
                } else {
                    0.0
                };
                let n50 = calc_n50(&lens);
                rows.push(StatRow {
                    file: p.clone(),
                    format: format!("{:?}", fmt).to_lowercase(),
                    records,
                    total_bases,
                    min_len,
                    max_len,
                    mean_len,
                    gc_pct,
                    n50,
                });
            }
            SeqFormat::Sam => {
                let content = std::fs::read_to_string(&p)?;
                let records = content.lines().filter(|l| !l.starts_with('@')).count();
                rows.push(StatRow {
                    file: p.clone(),
                    format: "sam".into(),
                    records,
                    total_bases: 0,
                    min_len: 0,
                    max_len: 0,
                    mean_len: 0.0,
                    gc_pct: 0.0,
                    n50: 0,
                });
            }
            SeqFormat::Bam | SeqFormat::Cram => {
                bail!("BAM/CRAM stats are not yet implemented in this build")
            }
        }
    }

    if args.json {
        println!("{}", serde_json::to_string_pretty(&rows)?);
    } else if args.tabular {
        println!("file\tformat\trecords\ttotal_bases\tmin_len\tmax_len\tmean_len\tgc_pct\tn50");
        for r in rows {
            println!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{}",
                r.file,
                r.format,
                r.records,
                r.total_bases,
                r.min_len,
                r.max_len,
                r.mean_len,
                r.gc_pct,
                r.n50
            );
        }
    } else {
        for r in rows {
            println!(
                "{} ({}) records={} bases={} min={} max={} mean={:.2} gc={:.2}% n50={}",
                r.file,
                r.format,
                r.records,
                r.total_bases,
                r.min_len,
                r.max_len,
                r.mean_len,
                r.gc_pct,
                r.n50
            );
        }
    }
    Ok(())
}

fn calc_n50(lens: &[usize]) -> usize {
    if lens.is_empty() {
        return 0;
    }
    let mut v = lens.to_vec();
    v.sort_unstable_by(|a, b| b.cmp(a));
    let sum: usize = v.iter().sum();
    let mut running = 0usize;
    for l in v {
        running += l;
        if running * 2 >= sum {
            return l;
        }
    }
    0
}
