use anyhow::{bail, Result};

use crate::formats::SeqRecord;

pub fn validate_pair_counts(r1: &[SeqRecord], r2: &[SeqRecord]) -> Result<()> {
    if r1.len() != r2.len() {
        bail!(
            "paired inputs have mismatched record counts: {} vs {}",
            r1.len(),
            r2.len()
        );
    }
    Ok(())
}

pub fn zip_pairs<'a>(
    r1: &'a [SeqRecord],
    r2: &'a [SeqRecord],
) -> Result<Vec<(&'a SeqRecord, &'a SeqRecord)>> {
    validate_pair_counts(r1, r2)?;
    Ok(r1.iter().zip(r2.iter()).collect())
}
