# seqit

`seqit` is a streaming-first Rust CLI toolkit for practical sequence-file manipulation in production bioinformatics pipelines.

## Current scope

This foundation release focuses on robust FASTA/FASTQ handling with paired-end correctness for pair-aware commands.

Implemented commands:

- `stats` – per-file and aggregate summary metrics.
- `seq` – general sequence filtering/transformation (length, quality, casing, reverse/comp, etc.).
- `fq2fa` – FASTQ to FASTA conversion.
- `grep` – pattern filtering on id/name/seq/qual, including pair-aware selection.
- `locate` – motif coordinate reporting.
- `sample` – random fixed-count or rate sampling.
- `rmdup` – duplicate removal/marking.
- `rename` – ID renaming (sequential, mapping-file, regex replacement).
- `sort` – sorting by ID/name/length/sequence.
- `shuffle` – deterministic random shuffling.
- `spike` – spike-in records from another dataset.
- `head` – keep first N or first proportion.
- `tail` – keep last N or last proportion.

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

- `-t, --threads` (set worker thread count; applies to all commands including `stats` and `spike`)
- `-o, --output`
- `-O, --output2`
- `-i, --input` (paired read 1 for pair-aware commands)
- `-I, --input2` (paired read 2 for pair-aware commands)
- `-s, --seed`
- `--format`
- `--compression`

## Examples

### stats

```bash
seqit stats reads.fq.gz
seqit stats reads_1.fq.gz reads_2.fq.gz                # pretty table with TOTAL row
seqit stats reads_1.fq.gz reads_2.fq.gz --per-file     # pretty table without TOTAL row
seqit stats reads.fq.gz -T                              # TSV output
seqit stats sample.fa --json
seqit stats sample.fa -t 8
```

`stats` output modes:

- default: human-readable pretty table with ASCII borders
- `-T, --tsv` (alias: `--tabular`): tab-separated output for scripting
- `--json`: structured JSON output

### seq

```bash
seqit seq reads.fa --min-len 75 --revcomp -o filtered.fa
```

Notes:
- `--min-qual/-q` and `--max-qual/-Q` are disabled by default with `-1` (not `-2`).
- Use `--name`, `--only-id`, or `--seq` for text-only projections.

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
seqit sample -i r1.fq -I r2.fq -r 0.1 -s 42 -o s1.fq -O s2.fq
```

### rmdup

```bash
seqit rmdup reads.fa --by seq --keep-first -o dedup.fa
seqit rmdup -i r1.fq -I r2.fq --by full --mark-dup -o d1.fq -O d2.fq
```

### rename

```bash
seqit rename reads.fa --prefix sample_ --start 1 --width 8 -o renamed.fa
seqit rename -i r1.fq -I r2.fq --prefix pair_ --keep-pair-suffix -o r1.new.fq -O r2.new.fq
seqit rename reads.fa -e 'lib_{n}_{prefix}' -p r -w 4 -o renamed.fa
seqit rename reads.fa --map-file id_map.tsv -o renamed.by_map.fa
seqit rename reads.fa --match-regex '^sample_(\d+)$' --replace 'S$1' -o renamed.by_regex.fa
```

`rename` modes:

- **Generate mode** (default): build IDs from `--prefix`, `--start`, `--width`, optional `--template`.
  - Template placeholders:
    - `{prefix}`: value from `--prefix`
    - `{n}`: zero-padded index based on `--start`/`--width`
- **Map mode**: `--map-file <TSV>`
  - File format: **exactly two tab-separated columns, no header**.
  - Column 1: existing ID, column 2: replacement ID.
  - Unmatched IDs are kept unchanged.
- **Regex mode**: `--match-regex <PATTERN> --replace <REPLACEMENT>`
  - Uses Rust regex syntax.
  - Replacement supports capture groups like `$1`, `$2`, etc.

Mode selection:
- `--mode auto` (default): picks map mode when `--map-file` is set, regex mode when `--match-regex/--replace` are set, otherwise generate mode.
- You can force with `--mode generate|map|regex`.

### sort

```bash
seqit sort reads.fa --by len --reverse -o sorted.fa
```

### shuffle

```bash
seqit shuffle reads.fa -s 42 -o shuffled.fa
seqit shuffle -i r1.fq -I r2.fq -s 42 -o shuf.r1.fq -O shuf.r2.fq
```

### spike

```bash
seqit spike -i target.fa -a inserts.fa -s 123 -o spiked.fa
seqit spike --input1 target.r1.fq -I target.r2.fq -a add.r1.fq -A add.r2.fq -s 123 -o out.r1.fq -O out.r2.fq
```

### head / tail

```bash
seqit head reads.fa -n 100 -o first100.fa
seqit tail reads.fq -p 0.1 -o last10pct.fq
seqit head -i r1.fq -I r2.fq -n 50000 -o h1.fq -O h2.fq
```

## Paired-end behavior

Pair-aware commands (`grep`, `sample`, `rmdup`, `rename`, `shuffle`, `spike`, `head`, `tail`) validate mate IDs and preserve pair lockstep.

- By default, invalid pairs fail fast with example read IDs and a total invalid count.
- Use `--allow-unpaired` to continue and skip invalid/unpaired records.
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
