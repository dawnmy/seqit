use anyhow::Result;
use num_format::{Locale, ToFormattedString};
use rayon::prelude::*;
use rust_htslib::bam::{self, Read};
use serde::Serialize;
use std::io::BufRead;

use crate::cli::StatsArgs;
use crate::formats::SeqFormat;
use crate::io;

#[derive(Clone, Copy)]
struct SeqAcc {
    records: usize,
    total_bases: usize,
    min_len: usize,
    max_len: usize,
    gc: usize,
}

impl SeqAcc {
    fn new() -> Self {
        Self {
            records: 0,
            total_bases: 0,
            min_len: usize::MAX,
            max_len: 0,
            gc: 0,
        }
    }

    fn with_record(len: usize, gc: usize) -> Self {
        Self {
            records: 1,
            total_bases: len,
            min_len: len,
            max_len: len,
            gc,
        }
    }

    fn merge(self, other: Self) -> Self {
        Self {
            records: self.records + other.records,
            total_bases: self.total_bases + other.total_bases,
            min_len: self.min_len.min(other.min_len),
            max_len: self.max_len.max(other.max_len),
            gc: self.gc + other.gc,
        }
    }
}

#[derive(Debug, Serialize, Clone)]
struct StatRow {
    file: String,
    format: String,
    seq_type: String,
    records: usize,
    total_bases: usize,
    min_len: usize,
    max_len: usize,
    mean_len: f64,
    gc_pct: f64,
    n50: usize,
    l50: usize,
    q1_len: usize,
    median_len: usize,
    q3_len: usize,
    n_pct: f64,
    q10_pct: Option<f64>,
    q20_pct: Option<f64>,
    q30_pct: Option<f64>,
}

#[derive(Clone, Copy)]
struct AllColumns {
    show_nx: bool,
    show_qx: bool,
}

#[derive(Default, Clone, Copy)]
struct SeqTypeCounts {
    alpha: usize,
    non_nuc_alpha: usize,
    u: usize,
    t: usize,
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
    let all_columns = choose_all_columns(&rows);

    if args.json {
        println!("{}", serde_json::to_string_pretty(&rows)?);
    } else if args.tsv {
        if args.all {
            print_tsv_all_header(all_columns);
        } else {
            println!(
                "file\tformat\ttype\trecords\ttotal_bases\tmin_len\tmax_len\tmean_len\tgc_pct\tn50"
            );
        }
        for r in rows {
            if args.all {
                println!("{}", tsv_all_row(&r, all_columns));
            } else {
                println!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{}",
                    r.file,
                    r.format,
                    r.seq_type,
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
    } else {
        if !args.per_file && rows.len() > 1 {
            rows.push(aggregate_row(&rows));
        }
        print_pretty_table(&rows, args.all, all_columns);
    }
    Ok(())
}

fn build_row(p: &str, args: &StatsArgs) -> Result<StatRow> {
    if p == "-" && matches!(args.format, crate::cli::FormatArg::Auto) {
        return build_seq_row_streaming(p, &args.format);
    }
    let fmt = SeqFormat::from_arg_or_detect(&args.format, Some(p))?;
    match fmt {
        SeqFormat::Fasta | SeqFormat::Fastq => build_seq_row_streaming(p, &args.format),
        SeqFormat::Sam => build_sam_row_streaming(p),
        SeqFormat::Bam | SeqFormat::Cram => build_hts_row(p, fmt),
    }
}

fn build_seq_row_streaming(file: &str, format_arg: &crate::cli::FormatArg) -> Result<StatRow> {
    let (fmt, br) = io::open_seq_reader(Some(file), format_arg, &crate::cli::CompressionArg::Auto)?;
    match fmt {
        SeqFormat::Fasta => build_fasta_row(file, br),
        SeqFormat::Fastq => build_fastq_row(file, br),
        _ => unreachable!("open_seq_reader only resolves fasta/fastq for stats sequence path"),
    }
}

fn build_fasta_row(file: &str, reader: impl BufRead) -> Result<StatRow> {
    let fa = bio::io::fasta::Reader::new(reader);
    let mut acc = SeqAcc::new();
    let mut lens = Vec::new();
    let mut n_bases = 0usize;
    let mut type_counts = SeqTypeCounts::default();

    for r in fa.records() {
        let r = r?;
        let seq = r.seq();
        let len = seq.len();
        acc = acc.merge(SeqAcc::with_record(
            len,
            bytecount::count(seq, b'G')
                + bytecount::count(seq, b'g')
                + bytecount::count(seq, b'C')
                + bytecount::count(seq, b'c'),
        ));
        n_bases += bytecount::count(seq, b'N') + bytecount::count(seq, b'n');
        update_seq_type_counts(seq, &mut type_counts);
        lens.push(len);
    }

    finalize_seq_row(
        file,
        SeqFormat::Fasta,
        acc,
        lens,
        n_bases,
        (None, None, None),
        type_counts,
    )
}

fn build_fastq_row(file: &str, reader: impl BufRead) -> Result<StatRow> {
    let fq = bio::io::fastq::Reader::new(reader);
    let mut acc = SeqAcc::new();
    let mut lens = Vec::new();
    let mut n_bases = 0usize;
    let mut q_ge_10 = 0usize;
    let mut q_ge_20 = 0usize;
    let mut q_ge_30 = 0usize;
    let mut qual_bases = 0usize;
    let mut type_counts = SeqTypeCounts::default();

    for r in fq.records() {
        let r = r?;
        let seq = r.seq();
        let len = seq.len();
        acc = acc.merge(SeqAcc::with_record(
            len,
            bytecount::count(seq, b'G')
                + bytecount::count(seq, b'g')
                + bytecount::count(seq, b'C')
                + bytecount::count(seq, b'c'),
        ));
        n_bases += bytecount::count(seq, b'N') + bytecount::count(seq, b'n');
        update_seq_type_counts(seq, &mut type_counts);
        lens.push(len);

        for q in r.qual().iter().copied() {
            let phred = q.saturating_sub(33);
            if phred >= 10 {
                q_ge_10 += 1;
            }
            if phred >= 20 {
                q_ge_20 += 1;
            }
            if phred >= 30 {
                q_ge_30 += 1;
            }
            qual_bases += 1;
        }
    }

    let q_pcts = (
        Some(calc_pct(q_ge_10, qual_bases)),
        Some(calc_pct(q_ge_20, qual_bases)),
        Some(calc_pct(q_ge_30, qual_bases)),
    );
    finalize_seq_row(
        file,
        SeqFormat::Fastq,
        acc,
        lens,
        n_bases,
        q_pcts,
        type_counts,
    )
}

fn finalize_seq_row(
    file: &str,
    fmt: SeqFormat,
    acc: SeqAcc,
    lens: Vec<usize>,
    n_bases: usize,
    q_pcts: (Option<f64>, Option<f64>, Option<f64>),
    type_counts: SeqTypeCounts,
) -> Result<StatRow> {
    let records = acc.records;
    let total_bases = acc.total_bases;
    let min_len = if records > 0 { acc.min_len } else { 0 };
    let max_len = acc.max_len;
    let mean_len = if records > 0 {
        total_bases as f64 / records as f64
    } else {
        0.0
    };
    let gc_pct = if total_bases > 0 {
        (acc.gc as f64) * 100.0 / (total_bases as f64)
    } else {
        0.0
    };
    let (n50, l50) = calc_nx_lx(&lens, 50.0);
    let (q1_len, median_len, q3_len) = calc_quartiles(&lens);
    let n_pct = calc_pct(n_bases, total_bases);
    let (q10_pct, q20_pct, q30_pct) = q_pcts;
    Ok(StatRow {
        file: file.to_string(),
        format: format!("{:?}", fmt).to_lowercase(),
        seq_type: infer_seq_type(type_counts),
        records,
        total_bases,
        min_len,
        max_len,
        mean_len,
        gc_pct,
        n50,
        l50,
        q1_len,
        median_len,
        q3_len,
        n_pct,
        q10_pct,
        q20_pct,
        q30_pct,
    })
}

fn build_sam_row_streaming(p: &str) -> Result<StatRow> {
    let r = io::open_reader(Some(p))?;
    let r = io::wrap_decompress(r, Some(p), &crate::cli::CompressionArg::Auto)?;
    let br = std::io::BufReader::new(r);
    let mut records = 0usize;
    for line in br.lines() {
        let line = line?;
        if !line.starts_with('@') {
            records += 1;
        }
    }
    Ok(StatRow {
        file: p.to_string(),
        format: "sam".into(),
        seq_type: "NA".into(),
        records,
        total_bases: 0,
        min_len: 0,
        max_len: 0,
        mean_len: 0.0,
        gc_pct: 0.0,
        n50: 0,
        l50: 0,
        q1_len: 0,
        median_len: 0,
        q3_len: 0,
        n_pct: 0.0,
        q10_pct: None,
        q20_pct: None,
        q30_pct: None,
    })
}

fn build_hts_row(p: &str, fmt: SeqFormat) -> Result<StatRow> {
    let mut reader = bam::Reader::from_path(p)?;
    let mut records = 0usize;
    let mut total_bases = 0usize;
    let mut min_len = usize::MAX;
    let mut max_len = 0usize;
    let mut gc = 0usize;
    let mut n_bases = 0usize;
    let mut q_ge_10 = 0usize;
    let mut q_ge_20 = 0usize;
    let mut q_ge_30 = 0usize;
    let mut qual_bases = 0usize;
    let mut lens = Vec::new();
    let mut type_counts = SeqTypeCounts::default();

    for rec in reader.records() {
        let rec = rec?;
        let seq = rec.seq().as_bytes();
        let len = seq.len();
        records += 1;
        total_bases += len;
        min_len = min_len.min(len);
        max_len = max_len.max(len);
        gc += bytecount::count(&seq, b'G')
            + bytecount::count(&seq, b'g')
            + bytecount::count(&seq, b'C')
            + bytecount::count(&seq, b'c');
        n_bases += bytecount::count(&seq, b'N') + bytecount::count(&seq, b'n');
        update_seq_type_counts(&seq, &mut type_counts);
        for q in rec.qual().iter() {
            let phred = q.saturating_sub(33);
            if phred >= 10 {
                q_ge_10 += 1;
            }
            if phred >= 20 {
                q_ge_20 += 1;
            }
            if phred >= 30 {
                q_ge_30 += 1;
            }
            qual_bases += 1;
        }
        lens.push(len);
    }

    let min_len = if records > 0 { min_len } else { 0 };
    let mean_len = if records > 0 {
        total_bases as f64 / records as f64
    } else {
        0.0
    };
    let gc_pct = if total_bases > 0 {
        (gc as f64) * 100.0 / (total_bases as f64)
    } else {
        0.0
    };

    let (n50, l50) = calc_nx_lx(&lens, 50.0);
    let (q1_len, median_len, q3_len) = calc_quartiles(&lens);
    Ok(StatRow {
        file: p.to_string(),
        format: format!("{:?}", fmt).to_lowercase(),
        seq_type: infer_seq_type(type_counts),
        records,
        total_bases,
        min_len,
        max_len,
        mean_len,
        gc_pct,
        n50,
        l50,
        q1_len,
        median_len,
        q3_len,
        n_pct: if total_bases > 0 {
            (n_bases as f64) * 100.0 / (total_bases as f64)
        } else {
            0.0
        },
        q10_pct: Some(calc_pct(q_ge_10, qual_bases)),
        q20_pct: Some(calc_pct(q_ge_20, qual_bases)),
        q30_pct: Some(calc_pct(q_ge_30, qual_bases)),
    })
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
    let l50 = rows.iter().map(|r| r.l50).max().unwrap_or(0);
    let q1_len = rows.iter().map(|r| r.q1_len).max().unwrap_or(0);
    let median_len = rows.iter().map(|r| r.median_len).max().unwrap_or(0);
    let q3_len = rows.iter().map(|r| r.q3_len).max().unwrap_or(0);
    let weighted_n_sum: f64 = rows.iter().map(|r| r.n_pct * r.total_bases as f64).sum();
    let n_pct = if total_bases > 0 {
        weighted_n_sum / total_bases as f64
    } else {
        0.0
    };
    let q10_pct = weighted_opt_pct(rows, |r| r.q10_pct);
    let q20_pct = weighted_opt_pct(rows, |r| r.q20_pct);
    let q30_pct = weighted_opt_pct(rows, |r| r.q30_pct);
    StatRow {
        file: "TOTAL".to_string(),
        format: "mixed".to_string(),
        seq_type: "mixed".to_string(),
        records,
        total_bases,
        min_len,
        max_len,
        mean_len,
        gc_pct,
        n50,
        l50,
        q1_len,
        median_len,
        q3_len,
        n_pct,
        q10_pct,
        q20_pct,
        q30_pct,
    }
}

fn print_pretty_table(rows: &[StatRow], all: bool, all_columns: AllColumns) {
    let headers = if all {
        let mut out = vec![
            "file".to_string(),
            "format".to_string(),
            "type".to_string(),
            "records".to_string(),
            "total_bases".to_string(),
            "min_len".to_string(),
            "max_len".to_string(),
            "mean_len".to_string(),
            "gc_pct".to_string(),
            "q1_len".to_string(),
            "median_len".to_string(),
            "q3_len".to_string(),
            "n_pct".to_string(),
        ];
        if all_columns.show_nx {
            out.extend(["n50".to_string(), "l50".to_string()]);
        }
        if all_columns.show_qx {
            out.extend([
                "q10_pct".to_string(),
                "q20_pct".to_string(),
                "q30_pct".to_string(),
            ]);
        }
        out
    } else {
        vec![
            "file".to_string(),
            "format".to_string(),
            "type".to_string(),
            "records".to_string(),
            "total_bases".to_string(),
            "min_len".to_string(),
            "max_len".to_string(),
            "mean_len".to_string(),
            "gc_pct".to_string(),
            "n50".to_string(),
        ]
    };
    let body: Vec<Vec<String>> = rows
        .iter()
        .map(|r| {
            if all {
                let mut out = vec![
                    r.file.clone(),
                    r.format.clone(),
                    r.seq_type.clone(),
                    r.records.to_formatted_string(&Locale::en),
                    r.total_bases.to_formatted_string(&Locale::en),
                    r.min_len.to_formatted_string(&Locale::en),
                    r.max_len.to_formatted_string(&Locale::en),
                    format!("{:.2}", r.mean_len),
                    format!("{:.2}%", r.gc_pct),
                    r.q1_len.to_formatted_string(&Locale::en),
                    r.median_len.to_formatted_string(&Locale::en),
                    r.q3_len.to_formatted_string(&Locale::en),
                    format!("{:.2}%", r.n_pct),
                ];
                if all_columns.show_nx {
                    out.extend([
                        display_usize(r.n50, !is_assembly(&r.format)),
                        display_usize(r.l50, !is_assembly(&r.format)),
                    ]);
                }
                if all_columns.show_qx {
                    out.extend([
                        display_opt_pct(r.q10_pct),
                        display_opt_pct(r.q20_pct),
                        display_opt_pct(r.q30_pct),
                    ]);
                }
                out
            } else {
                vec![
                    r.file.clone(),
                    r.format.clone(),
                    r.seq_type.clone(),
                    r.records.to_formatted_string(&Locale::en),
                    r.total_bases.to_formatted_string(&Locale::en),
                    r.min_len.to_formatted_string(&Locale::en),
                    r.max_len.to_formatted_string(&Locale::en),
                    format!("{:.2}", r.mean_len),
                    format!("{:.2}%", r.gc_pct),
                    r.n50.to_formatted_string(&Locale::en),
                ]
            }
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

fn choose_all_columns(rows: &[StatRow]) -> AllColumns {
    let has_assembly = rows.iter().any(|r| is_assembly(&r.format));
    let has_reads = rows.iter().any(|r| is_read_or_alignment(&r.format));
    match (has_assembly, has_reads) {
        (true, true) => AllColumns {
            show_nx: true,
            show_qx: true,
        },
        (true, false) => AllColumns {
            show_nx: true,
            show_qx: false,
        },
        (false, true) => AllColumns {
            show_nx: false,
            show_qx: true,
        },
        (false, false) => AllColumns {
            show_nx: false,
            show_qx: false,
        },
    }
}

fn is_assembly(format: &str) -> bool {
    format == "fasta"
}

fn is_read_or_alignment(format: &str) -> bool {
    matches!(format, "fastq" | "sam" | "bam" | "cram")
}

fn print_tsv_all_header(columns: AllColumns) {
    let mut headers = vec![
        "file",
        "format",
        "type",
        "records",
        "total_bases",
        "min_len",
        "max_len",
        "mean_len",
        "gc_pct",
        "q1_len",
        "median_len",
        "q3_len",
        "n_pct",
    ];
    if columns.show_nx {
        headers.extend(["n50", "l50"]);
    }
    if columns.show_qx {
        headers.extend(["q10_pct", "q20_pct", "q30_pct"]);
    }
    println!("{}", headers.join("\t"));
}

fn tsv_all_row(r: &StatRow, columns: AllColumns) -> String {
    let mut fields = vec![
        r.file.clone(),
        r.format.clone(),
        r.seq_type.clone(),
        r.records.to_string(),
        r.total_bases.to_string(),
        r.min_len.to_string(),
        r.max_len.to_string(),
        format!("{:.3}", r.mean_len),
        format!("{:.3}", r.gc_pct),
        r.q1_len.to_string(),
        r.median_len.to_string(),
        r.q3_len.to_string(),
        format!("{:.3}", r.n_pct),
    ];
    if columns.show_nx {
        fields.extend([
            display_usize(r.n50, !is_assembly(&r.format)),
            display_usize(r.l50, !is_assembly(&r.format)),
        ]);
    }
    if columns.show_qx {
        fields.extend([
            display_opt_pct(r.q10_pct),
            display_opt_pct(r.q20_pct),
            display_opt_pct(r.q30_pct),
        ]);
    }
    fields.join("\t")
}

fn update_seq_type_counts(seq: &[u8], counts: &mut SeqTypeCounts) {
    for b in seq.iter().copied() {
        if !b.is_ascii_alphabetic() {
            continue;
        }
        counts.alpha += 1;
        let up = b.to_ascii_uppercase();
        match up {
            b'U' => counts.u += 1,
            b'T' => counts.t += 1,
            _ => {}
        }
        if !matches!(up, b'A' | b'C' | b'G' | b'T' | b'U' | b'N') {
            counts.non_nuc_alpha += 1;
        }
    }
}

fn infer_seq_type(counts: SeqTypeCounts) -> String {
    if counts.alpha == 0 {
        return "NA".to_string();
    }
    let non_nuc_ratio = counts.non_nuc_alpha as f64 / counts.alpha as f64;
    if non_nuc_ratio > 0.05 {
        return "protein".to_string();
    }
    if counts.u > 0 && counts.t == 0 {
        return "RNA".to_string();
    }
    "DNA".to_string()
}

fn calc_nx_lx(lens: &[usize], x: f64) -> (usize, usize) {
    if lens.is_empty() {
        return (0, 0);
    }
    let mut v = lens.to_vec();
    v.sort_unstable_by(|a, b| b.cmp(a));
    let sum: usize = v.iter().sum();
    let threshold = (sum as f64) * (x / 100.0);
    let mut running = 0usize;
    for (idx, l) in v.iter().copied().enumerate() {
        running += l;
        if (running as f64) >= threshold {
            return (l, idx + 1);
        }
    }
    (0, 0)
}

fn calc_quartiles(lens: &[usize]) -> (usize, usize, usize) {
    if lens.is_empty() {
        return (0, 0, 0);
    }
    let mut v = lens.to_vec();
    v.sort_unstable();
    (
        percentile(&v, 25.0),
        percentile(&v, 50.0),
        percentile(&v, 75.0),
    )
}

fn percentile(sorted: &[usize], pct: f64) -> usize {
    if sorted.is_empty() {
        return 0;
    }
    let n = sorted.len();
    let idx = (((pct / 100.0) * ((n - 1) as f64)).round()) as usize;
    sorted[idx]
}

fn calc_pct(numerator: usize, denominator: usize) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        (numerator as f64) * 100.0 / (denominator as f64)
    }
}

fn weighted_opt_pct(rows: &[StatRow], pick: fn(&StatRow) -> Option<f64>) -> Option<f64> {
    let mut weighted = 0.0f64;
    let mut denom = 0usize;
    for r in rows {
        if let Some(v) = pick(r) {
            weighted += v * (r.total_bases as f64);
            denom += r.total_bases;
        }
    }
    if denom == 0 {
        None
    } else {
        Some(weighted / denom as f64)
    }
}

fn display_usize(value: usize, hide: bool) -> String {
    if hide {
        "NA".to_string()
    } else {
        value.to_formatted_string(&Locale::en)
    }
}

fn display_opt_pct(value: Option<f64>) -> String {
    value
        .map(|v| format!("{v:.2}%"))
        .unwrap_or_else(|| "NA".to_string())
}
