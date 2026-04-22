use anyhow::{bail, Result};
use num_format::{Locale, ToFormattedString};
use rayon::prelude::*;
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
    let mut rows = inputs
        .par_iter()
        .map(|p| build_row(p, &args))
        .collect::<Result<Vec<_>>>()?;
    rows.sort_unstable_by(|a, b| a.file.cmp(&b.file));

    if args.json {
        println!("{}", serde_json::to_string_pretty(&rows)?);
    } else if args.tsv {
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
        if !args.per_file && rows.len() > 1 {
            rows.push(aggregate_row(&rows));
        }
        print_pretty_table(&rows);
    }
    Ok(())
}

fn build_row(p: &str, args: &StatsArgs) -> Result<StatRow> {
    let fmt = SeqFormat::from_arg(&args.format).unwrap_or(SeqFormat::detect(Some(p))?);
    match fmt {
        SeqFormat::Fasta | SeqFormat::Fastq => {
            let recs = io::read_records(Some(p), fmt, &crate::cli::CompressionArg::Auto)?;
            let lens: Vec<usize> = recs.par_iter().map(|r| r.seq.len()).collect();
            let records = lens.len();
            let total_bases: usize = lens.par_iter().sum();
            let min_len = *lens.par_iter().min().unwrap_or(&0);
            let max_len = *lens.par_iter().max().unwrap_or(&0);
            let mean_len = if records > 0 {
                total_bases as f64 / records as f64
            } else {
                0.0
            };
            let gc = recs
                .par_iter()
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
            Ok(StatRow {
                file: p.to_string(),
                format: format!("{:?}", fmt).to_lowercase(),
                records,
                total_bases,
                min_len,
                max_len,
                mean_len,
                gc_pct,
                n50,
            })
        }
        SeqFormat::Sam => {
            let content = std::fs::read_to_string(p)?;
            let records = content.lines().filter(|l| !l.starts_with('@')).count();
            Ok(StatRow {
                file: p.to_string(),
                format: "sam".into(),
                records,
                total_bases: 0,
                min_len: 0,
                max_len: 0,
                mean_len: 0.0,
                gc_pct: 0.0,
                n50: 0,
            })
        }
        SeqFormat::Bam | SeqFormat::Cram => {
            bail!("BAM/CRAM stats are not yet implemented in this build")
        }
    }
}

fn aggregate_row(rows: &[StatRow]) -> StatRow {
    let records: usize = rows.iter().map(|r| r.records).sum();
    let total_bases: usize = rows.iter().map(|r| r.total_bases).sum();
    let min_len = rows.iter().map(|r| r.min_len).min().unwrap_or(0);
    let max_len = rows.iter().map(|r| r.max_len).max().unwrap_or(0);
    let mean_len = if records > 0 {
        total_bases as f64 / records as f64
    } else {
        0.0
    };
    let weighted_gc_sum: f64 = rows.iter().map(|r| r.gc_pct * r.total_bases as f64).sum();
    let gc_pct = if total_bases > 0 {
        weighted_gc_sum / total_bases as f64
    } else {
        0.0
    };
    let n50 = rows.iter().map(|r| r.n50).max().unwrap_or(0);
    StatRow {
        file: "TOTAL".to_string(),
        format: "mixed".to_string(),
        records,
        total_bases,
        min_len,
        max_len,
        mean_len,
        gc_pct,
        n50,
    }
}

fn print_pretty_table(rows: &[StatRow]) {
    let headers = vec![
        "file".to_string(),
        "format".to_string(),
        "records".to_string(),
        "total_bases".to_string(),
        "min_len".to_string(),
        "max_len".to_string(),
        "mean_len".to_string(),
        "gc_pct".to_string(),
        "n50".to_string(),
    ];
    let body: Vec<Vec<String>> = rows
        .iter()
        .map(|r| {
            vec![
                r.file.clone(),
                r.format.clone(),
                r.records.to_formatted_string(&Locale::en),
                r.total_bases.to_formatted_string(&Locale::en),
                r.min_len.to_formatted_string(&Locale::en),
                r.max_len.to_formatted_string(&Locale::en),
                format!("{:.2}", r.mean_len),
                format!("{:.2}%", r.gc_pct),
                r.n50.to_formatted_string(&Locale::en),
            ]
        })
        .collect();
    let mut widths: Vec<usize> = headers.iter().map(|h| h.len()).collect();
    for row in &body {
        for (i, cell) in row.iter().enumerate() {
            widths[i] = widths[i].max(cell.len());
        }
    }
    let border = format!(
        "+-{}-+",
        widths
            .iter()
            .map(|w| "-".repeat(*w))
            .collect::<Vec<_>>()
            .join("-+-")
    );
    println!("{border}");
    println!(
        "| {} |",
        headers
            .iter()
            .enumerate()
            .map(|(i, h)| format!("{h:<width$}", width = widths[i]))
            .collect::<Vec<_>>()
            .join(" | ")
    );
    println!("{border}");
    for row in body {
        println!(
            "| {} |",
            row.iter()
                .enumerate()
                .map(|(i, c)| format!("{c:<width$}", width = widths[i]))
                .collect::<Vec<_>>()
                .join(" | ")
        );
    }
    println!("{border}");
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
