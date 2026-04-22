# seqit

`seqit` is a streaming-first Rust CLI toolkit for practical sequence-file manipulation in production bioinformatics pipelines.

## Current scope

This foundation release focuses on robust FASTA/FASTQ handling with paired-end correctness for pair-aware commands.

Implemented commands:

- `stats`
- `seq`
- `fq2fa`
- `grep`
- `locate`
- `sample`
- `rmdup`
- `rename`
- `sort`
- `shuffle`
- `spike`

SAM parsing for `stats` (basic count) is available, while BAM/CRAM support is scaffolded for upcoming work.

## Installation

```bash
cargo build --release
./target/release/seqit --help
```

## Supported formats and compression

Input/output formats:

- FASTA
- FASTQ
- SAM (stats only, text)
- BAM/CRAM (detected, but not yet implemented in commands)

Compression:

- none
- gzip (`.gz`)
- xz (`.xz`)

## CLI design conventions

Common options are unified across commands where relevant:

- `-t, --threads`
- `-o, --output`
- `-O, --output2`
- `-1, --in1`
- `-2, --in2`
- `-s, --seed`
- `--format`
- `--compression`

## Examples

### stats

```bash
seqit stats reads.fq.gz --tabular
seqit stats sample.fa --json
```

### seq

```bash
seqit seq reads.fa --min-len 75 --revcomp -o filtered.fa
```

### fq2fa

```bash
seqit fq2fa reads.fq.gz -o reads.fa.gz --compression gz
```

### grep

```bash
seqit grep reads.fq -p ACTG --by seq -o matched.fq
seqit grep --in1 r1.fq --in2 r2.fq -p adapter --pair-mode any -o out.r1.fq -O out.r2.fq
```

### locate

```bash
seqit locate reads.fa -p ACGT --all
seqit locate reads.fa -p 'A[CT]G' --regex --bed
```

### sample

```bash
seqit sample reads.fq -n 100000 -s 42 -o sub.fq
seqit sample --in1 r1.fq --in2 r2.fq -r 0.1 -s 42 -o s1.fq -O s2.fq
```

### rmdup

```bash
seqit rmdup reads.fa --by seq --keep-first -o dedup.fa
seqit rmdup --in1 r1.fq --in2 r2.fq --by full --mark-dup -o d1.fq -O d2.fq
```

### rename

```bash
seqit rename reads.fa --prefix sample_ --start 1 --width 8 -o renamed.fa
seqit rename --in1 r1.fq --in2 r2.fq --prefix pair_ --keep-pair-suffix -o r1.new.fq -O r2.new.fq
```

### sort

```bash
seqit sort reads.fa --by len --reverse -o sorted.fa
```

### shuffle

```bash
seqit shuffle reads.fa -s 42 -o shuffled.fa
seqit shuffle --in1 r1.fq --in2 r2.fq -s 42 -o shuf.r1.fq -O shuf.r2.fq
```

### spike

```bash
seqit spike -i target.fa -a inserts.fa -s 123 -o spiked.fa
seqit spike -1 target.r1.fq -2 target.r2.fq -a add.r1.fq -A add.r2.fq -s 123 -o out.r1.fq -O out.r2.fq
```

## Paired-end behavior

Pair-aware commands (`grep`, `sample`, `rmdup`, `rename`, `shuffle`, `spike`) validate mate counts and preserve pair lockstep.

- Any mismatch in paired input counts fails fast.
- Pair operations are done at pair granularity (not independent mates).

## Design notes

- Streaming-oriented I/O abstraction with format and compression auto-detection.
- Clean command module separation (`src/commands/*`) to simplify extension.
- Deterministic RNG for reproducible sampling/shuffling/spike-in.

## Tests

Run:

```bash
cargo test
```
