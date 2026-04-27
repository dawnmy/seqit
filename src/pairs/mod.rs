use anyhow::{bail, Result};

use crate::formats::SeqRecord;

pub fn prepare_paired_records(
    mut r1: Vec<SeqRecord>,
    mut r2: Vec<SeqRecord>,
    allow_unpaired: bool,
) -> Result<(Vec<SeqRecord>, Vec<SeqRecord>)> {
    let mut invalid_preview: Vec<String> = Vec::new();
    let mut invalid_count = 0usize;
    let min_len = r1.len().min(r2.len());
    for i in 0..min_len {
        let a = pair_key(&r1[i].id);
        let b = pair_key(&r2[i].id);
        if a != b {
            invalid_count += 1;
            if invalid_preview.len() < 10 {
                invalid_preview.push(format!("pair {}: R1={} R2={}", i + 1, r1[i].id, r2[i].id));
            }
        }
    }
    if r1.len() > min_len {
        invalid_count += r1.len() - min_len;
        for rec in r1.iter().skip(min_len).take(10 - invalid_preview.len()) {
            invalid_preview.push(format!("unmatched in R1: {}", rec.id));
        }
    }
    if r2.len() > min_len {
        invalid_count += r2.len() - min_len;
        for rec in r2.iter().skip(min_len).take(10 - invalid_preview.len()) {
            invalid_preview.push(format!("unmatched in R2: {}", rec.id));
        }
    }
    if invalid_count > 0 {
        let preview = invalid_preview.join("; ");
        if !allow_unpaired {
            bail!(
                "paired inputs contain {invalid_count} invalid/unpaired records; examples: {preview}. rerun with --allow-unpaired to continue"
            );
        }
        eprintln!(
            "warning: paired inputs contain {invalid_count} invalid/unpaired records; examples: {preview}. unpaired records will be skipped"
        );
    }

    let mut out1 = Vec::new();
    let mut out2 = Vec::new();
    for (a, b) in r1.drain(..min_len).zip(r2.drain(..min_len)) {
        if pair_key(&a.id) == pair_key(&b.id) {
            out1.push(a);
            out2.push(b);
        }
    }
    Ok((out1, out2))
}

fn pair_key(id: &str) -> &str {
    id.strip_suffix("/1")
        .or_else(|| id.strip_suffix("/2"))
        .unwrap_or(id)
}
