# seqit

`seqit` is a streaming-first Rust CLI toolkit for practical FASTA/FASTQ sequence-file operations in production pipelines.

## Installation

Download the latest precompiled binary for your platform from the [GitHub releases page](https://github.com/dawnmy/seqit/releases), or build from source as described below.

```bash
cargo build --release
./target/release/seqit --help
```

## Quick start

```bash
# Basic stats
seqit stats reads.fq.gz

# Filter by length and reverse-complement
seqit seq reads.fa --min-len 100 --revcomp -o filtered.fa

# Convert FASTQ to FASTA
seqit fq2fa reads.fq.gz -o reads.fa.gz --compression gz
```

## Input semantics (single-end vs paired-end)

For pair-aware commands (`grep`, `sample`, `rmdup`, `rename`, `shuffle`, `spike`, `head`, `tail`):

- **Single-end mode**: provide one positional `INPUT` (or stream from stdin when omitted).
- **Paired-end mode**: provide **both mates together** (`-i/-1/--in1` style + `-I/-2/--in2` style, depending on command).
- Mixed/partial mate flags are invalid.
- Use `-o/--output` and `-O/--output2` for paired outputs.
- `--allow-unpaired` keeps processing and skips invalid/unpaired records.

---

## Common options

Many commands share these options:

- `INPUT` (positional): single-end/single-file input.
- `-o, --output <FILE>`: primary output (`-` means stdout).
- `-O, --output2 <FILE>`: mate-2 output in paired mode.
- `--format <auto|fasta|fastq|sam|bam|cram>`: IO format override.
- `--compression <auto|none|gz|xz>`: output compression.
- `-t, --threads <N>`: worker thread count.
- `--quiet`: suppress non-error logs.

## Formats and compression

Supported formats:

- FASTA
- FASTQ
- SAM (`stats` supports basic counting)
- BAM/CRAM are scaffolded/detected for future support

Compression codecs:

- none
- gzip (`.gz`)
- xz (`.xz`)

---

## Command reference

## `seqit stats`
Compute summary statistics for one or more files.

### Key options
- `INPUTS...`: one or more input files.
- `-T, --tsv` / `--tabular`: TSV output.
- `-j, --json`: JSON output.
- `-a, --all`: extended metrics.
- `-p, --per-file`: suppress TOTAL row.
- `--format`, `--threads`.

### Examples
```bash
seqit stats reads.fq.gz
seqit stats r1.fq.gz r2.fq.gz --per-file
seqit stats assembly.fa -a
seqit stats sample.fa --json
```

## `seqit seq`
Filter/transform sequences and reads.

### Key options
- Length filters: `--min-len`, `--max-len`.
- Orientation: `--rev`, `--comp`, `--revcomp`.
- Case/gaps: `--upper`, `--lower`, `--remove-gaps`, `--gap-letters`.
- Quality filters: `--min-qual`, `--max-qual`, `--qual-ascii-base`.
- Projection/validation: `--name`, `--full-name`, `--only-id`, `--seq`, `--validate-seq`.
- Display: `--color`.

### Examples
```bash
seqit seq reads.fa --min-len 100 --max-len 2000 -o filtered.fa
seqit seq reads.fq --min-qual 20 --revcomp -o cleaned.fq
seqit seq reads.fa --name --full-name
```

## `seqit fq2fa`
Convert FASTQ to FASTA.

### Examples
```bash
seqit fq2fa reads.fq -o reads.fa
seqit fq2fa reads.fq.gz -o reads.fa.gz --compression gz
```

## `seqit grep`
Search/filter by ID, name, sequence, or quality.

### Key options
- Patterns: `-p/--pattern`, `-f/--pattern-file`.
- Search field: `-b/--by <id|name|seq|qual>`.
- Matching mode: `--regex`, `--ignore-case`, `--invert`.
- Output behavior: `--count`, `--only-names`.
- Paired mode: `--in1`, `--in2`, `--pair-mode <any|both>`, `--output2`, `--allow-unpaired`.

### Examples
```bash
seqit grep reads.fq -p ACTG --by seq -o matched.fq
seqit grep reads.fq --pattern-file ids.txt --by id --invert
seqit grep --in1 r1.fq --in2 r2.fq -p adapter --pair-mode both -o out.r1.fq -O out.r2.fq
```

## `seqit locate`
Locate motifs and report coordinates.

### Key options
- `-p/--pattern`, `-f/--pattern-file`
- `--regex`, `--ignore-case`
- `-a/--all`: report all matches per record.
- `-B/--bed`: BED-like output.

### Examples
```bash
seqit locate reads.fa -p ACGT
seqit locate reads.fa -p 'A[CT]G' --regex --all --bed
seqit locate reads.fa --pattern-file motifs.txt --ignore-case
```

## `seqit sample`
Random sampling by fixed count or rate.

### Key options
- Single-end positional `INPUT`.
- Paired flags: `-i/--input|--in1`, `-I/--input2|--in2`.
- `-n/--num`, `-r/--rate`.
- `-s/--seed`.
- `--progress`, `--allow-unpaired`.

### Examples
```bash
seqit sample reads.fq -n 100000 -s 42 -o sub.fq
seqit sample reads.fq -r 0.1 -s 42 -o sub.fq
seqit sample -i r1.fq -I r2.fq -n 50000 -o sub.r1.fq -O sub.r2.fq
```

## `seqit rmdup`
Remove/mark duplicates (single or paired).

### Key options
- Duplicate key: `-b/--by <id|seq|full>`.
- Retention: `--keep-first` or `--keep-last`.
- Annotation: `--count-dup`, `--mark-dup`.
- Paired options: `-i/-I`, `--in1/--in2`, `-O`, `--allow-unpaired`.

### Examples
```bash
seqit rmdup reads.fa --by seq --keep-first -o dedup.fa
seqit rmdup reads.fa --by full --keep-last -o dedup.fa
seqit rmdup -i r1.fq -I r2.fq --mark-dup -o d1.fq -O d2.fq
```

## `seqit rename`
Rename IDs with generation, mapping, or regex replacement.

### Key options
- Generate mode: `--prefix`, `--start`, `--width`, `--template`.
- Map mode: `--map-file <TSV>` (`old_id<TAB>new_id`, no header).
- Regex mode: `--match-regex`, `--replace`.
- Mode control: `--mode <auto|generate|map|regex>`.
- Paired behavior: `--keep-pair-suffix`, paired mate flags, `--allow-unpaired`.

### Template placeholders
- `{prefix}` → value of `--prefix`
- `{n}` → zero-padded index from `--start`/`--width`

### Examples
```bash
seqit rename reads.fa --prefix sample_ --start 1 --width 8 -o renamed.fa
seqit rename reads.fa --template 'lib_{prefix}_{n}' --prefix r --width 4 -o renamed.fa
seqit rename reads.fa --map-file id_map.tsv -o renamed.by_map.fa
seqit rename reads.fa --match-regex '^sample_(\d+)$' --replace 'S$1' -o renamed.by_regex.fa
seqit rename -i r1.fq -I r2.fq --prefix pair_ --keep-pair-suffix -o r1.new.fq -O r2.new.fq
```

## `seqit sort`
Sort records by key.

### Key options
- `-b/--by <id|name|len|seq>`.
- `-n/--numeric`.
- `-r/--reverse`.
- `-m/--mem` (spill threshold), `-T/--tmp-dir`.

### Examples
```bash
seqit sort reads.fa --by id -o sorted.fa
seqit sort reads.fa --by len --numeric --reverse -o sorted.fa
seqit sort reads.fa --by seq --mem 512M --tmp-dir /tmp -o sorted.fa
```

## `seqit shuffle`
Deterministic random shuffling.

### Key options
- `-s/--seed`.
- `-m/--mem`, `-T/--tmp-dir` for spill behavior.
- Paired mate flags and `--allow-unpaired`.

### Examples
```bash
seqit shuffle reads.fa -s 42 -o shuffled.fa
seqit shuffle reads.fq --mem 1G --tmp-dir /tmp -o shuffled.fq
seqit shuffle -i r1.fq -I r2.fq -s 42 -o shuf.r1.fq -O shuf.r2.fq
```

## `seqit spike`
Spike records from an additional dataset into a main dataset.

### Key options
- Main input: positional `INPUT` (single-end) or `-i/-I` mates (paired).
- Spike input: `-a/--add` (required), plus `-A/--add2` for paired.
- `-s/--seed`, `-o/-O`, `--format`, `--compression`, `--threads`, `--allow-unpaired`.

### Examples
```bash
seqit spike target.fa -a inserts.fa -s 123 -o spiked.fa
seqit spike target.fq -a add.fq -o spiked.fq
seqit spike -i target.r1.fq -I target.r2.fq -a add.r1.fq -A add.r2.fq -s 123 -o out.r1.fq -O out.r2.fq
```

## `seqit head`
Keep first N records/pairs or first proportion.

### Key options
- `-n/--num`, `-p/--proportion`
- Single-end positional `INPUT` or paired `-i/-I`
- `-O`, `--allow-unpaired`

### Examples
```bash
seqit head reads.fa -n 100 -o first100.fa
seqit head reads.fq -p 0.1 -o first10pct.fq
seqit head -i r1.fq -I r2.fq -n 50000 -o h1.fq -O h2.fq
```

## `seqit tail`
Keep last N records/pairs or last proportion.

### Key options
- `-n/--num`, `-p/--proportion`
- Single-end positional `INPUT` or paired `-i/-I`
- `-O`, `--allow-unpaired`

### Examples
```bash
seqit tail reads.fa -n 100 -o last100.fa
seqit tail reads.fq -p 0.1 -o last10pct.fq
seqit tail -i r1.fq -I r2.fq -n 50000 -o t1.fq -O t2.fq
```

---

## Testing

```bash
cargo test
```

### Acknowledgment
The FASTX parser used in `seqit` is based on the implementation from [shenwei356/fastx](https://github.com/shenwei356/fastx).
