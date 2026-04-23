use clap::{ArgAction, Parser, Subcommand, ValueEnum};

#[derive(Debug, Clone, ValueEnum)]
pub enum FormatArg {
    Auto,
    Fasta,
    Fastq,
    Sam,
    Bam,
    Cram,
}

#[derive(Debug, Clone, ValueEnum)]
pub enum CompressionArg {
    Auto,
    None,
    Gz,
    Xz,
}

#[derive(Debug, Parser)]
#[command(
    name = "seqit",
    version,
    about = "High-performance streaming toolkit for sequence files"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    #[command(about = "Compute summary statistics for sequence files")]
    Stats(StatsArgs),
    #[command(about = "Filter and transform sequences")]
    Seq(SeqArgs),
    #[command(about = "Convert FASTQ to FASTA")]
    Fq2fa(Fq2faArgs),
    #[command(about = "Search records by id/name/sequence/quality")]
    Grep(GrepArgs),
    #[command(about = "Locate sequence motifs and report coordinates")]
    Locate(LocateArgs),
    #[command(about = "Randomly sample records or pairs")]
    Sample(SampleArgs),
    #[command(about = "Remove duplicate records or read pairs")]
    Rmdup(RmdupArgs),
    #[command(about = "Rename record identifiers")]
    Rename(RenameArgs),
    #[command(about = "Sort records by id/name/length/sequence")]
    Sort(SortArgs),
    #[command(about = "Shuffle records or read pairs deterministically")]
    Shuffle(ShuffleArgs),
    #[command(about = "Spike reads/sequences from one dataset into another")]
    Spike(SpikeArgs),
    #[command(about = "Print the first N records/pairs or first proportion")]
    Head(HeadArgs),
    #[command(about = "Print the last N records/pairs or last proportion")]
    Tail(TailArgs),
}

impl Commands {
    pub fn threads(&self) -> Option<usize> {
        match self {
            Commands::Stats(a) => a.threads,
            Commands::Seq(a) => a.io.threads,
            Commands::Fq2fa(a) => a.io.threads,
            Commands::Grep(a) => a.io.threads,
            Commands::Locate(a) => a.io.threads,
            Commands::Sample(a) => a.io.threads,
            Commands::Rmdup(a) => a.io.threads,
            Commands::Rename(a) => a.io.threads,
            Commands::Sort(a) => a.io.threads,
            Commands::Shuffle(a) => a.io.threads,
            Commands::Spike(a) => a.threads,
            Commands::Head(a) => a.io.threads,
            Commands::Tail(a) => a.io.threads,
        }
    }
}

#[derive(Debug, clap::Args, Clone)]
pub struct CommonIoArgs {
    #[arg(
        value_name = "INPUT",
        help = "Positional single-end/single-file input (default: stdin). Pair-aware commands enter paired mode only when both mate flags are provided."
    )]
    pub input: Option<String>,
    #[arg(
        short = 'o',
        long = "output",
        default_value = "-",
        help = "Output file ('-' for stdout)"
    )]
    pub output: String,
    #[arg(
        long = "format",
        value_enum,
        default_value = "auto",
        help = "Input/output format"
    )]
    pub format: FormatArg,
    #[arg(
        long = "compression",
        value_enum,
        default_value = "auto",
        help = "Output compression codec"
    )]
    pub compression: CompressionArg,
    #[arg(
        short = 't',
        long = "threads",
        help = "Number of worker threads to use"
    )]
    pub threads: Option<usize>,
    #[arg(long = "quiet", action = ArgAction::SetTrue, help = "Suppress non-error logging")]
    pub quiet: bool,
}

#[derive(Debug, clap::Args)]
#[command(
    after_help = "Examples:\n  seqit stats reads.fq.gz\n  seqit stats r1.fq.gz r2.fq.gz --per-file\n  seqit stats assembly.fa -a --json"
)]
pub struct StatsArgs {
    #[arg(value_name = "INPUTS", help = "One or more input files")]
    pub inputs: Vec<String>,
    #[arg(
        short = 'T',
        long = "tsv",
        visible_alias = "tabular",
        action = ArgAction::SetTrue,
        help = "Output machine-friendly TSV instead of pretty table"
    )]
    pub tsv: bool,
    #[arg(short = 'j', long = "json", action = ArgAction::SetTrue, help = "Output JSON")]
    pub json: bool,
    #[arg(
        short = 'a',
        long = "all",
        action = ArgAction::SetTrue,
        help = "Enable extended assembly/read metrics (quartiles, L50, N%, and FASTQ Q10/Q20/Q30)"
    )]
    pub all: bool,
    #[arg(
        short = 'p',
        long = "per-file",
        action = ArgAction::SetTrue,
        help = "When multiple files are given, do not append a TOTAL summary row"
    )]
    pub per_file: bool,
    #[arg(
        long = "format",
        value_enum,
        default_value = "auto",
        help = "Input format"
    )]
    pub format: FormatArg,
    #[arg(
        short = 't',
        long = "threads",
        help = "Number of worker threads to use"
    )]
    pub threads: Option<usize>,
}

#[derive(Debug, clap::Args)]
#[command(
    after_help = "Examples:\n  seqit seq reads.fa --min-len 100 --max-len 500 -o filtered.fa\n  seqit seq reads.fq --min-qual 20 --revcomp -o cleaned.fq\n  seqit seq reads.fa --name --full-name"
)]
pub struct SeqArgs {
    #[command(flatten)]
    pub io: CommonIoArgs,
    #[arg(
        short = 'm',
        long = "min-len",
        help = "Keep only records with length >= N"
    )]
    pub min_len: Option<usize>,
    #[arg(
        short = 'M',
        long = "max-len",
        help = "Keep only records with length <= N"
    )]
    pub max_len: Option<usize>,
    #[arg(short = 'r', long = "rev", action = ArgAction::SetTrue, help = "Reverse sequence and quality order")]
    pub rev: bool,
    #[arg(short = 'c', long = "comp", action = ArgAction::SetTrue, help = "Complement sequence bases")]
    pub comp: bool,
    #[arg(short = 'R', long = "revcomp", action = ArgAction::SetTrue, help = "Reverse-complement sequence and reverse quality")]
    pub revcomp: bool,
    #[arg(short = 'u', long = "upper", action = ArgAction::SetTrue, help = "Convert sequence to uppercase")]
    pub upper: bool,
    #[arg(short = 'l', long = "lower", action = ArgAction::SetTrue, help = "Convert sequence to lowercase")]
    pub lower: bool,
    #[arg(
        short = 'k',
        long = "color",
        action = ArgAction::SetTrue,
        help = "Colorize sequence output with ANSI codes"
    )]
    pub color: bool,
    #[arg(
        short = 'G',
        long = "gap-letters",
        default_value = "- \t.",
        help = "Gap letters to remove with -g/--remove-gaps"
    )]
    pub gap_letters: String,
    #[arg(
        short = 'Q',
        long = "max-qual",
        default_value_t = -1.0,
        help = "Only keep reads with average quality <= N; -1 disables"
    )]
    pub max_qual: f64,
    #[arg(
        short = 'q',
        long = "min-qual",
        default_value_t = -1.0,
        help = "Only keep reads with average quality >= N; -1 disables"
    )]
    pub min_qual: f64,
    #[arg(
        short = 'n',
        long = "name",
        action = ArgAction::SetTrue,
        help = "Only print sequence IDs"
    )]
    pub name: bool,
    #[arg(
        long = "full-name",
        action = ArgAction::SetTrue,
        help = "With -n/--name, print full headers (ID + description) instead of IDs only"
    )]
    pub full_name: bool,
    #[arg(
        short = 'i',
        long = "only-id",
        action = ArgAction::SetTrue,
        help = "Only print IDs (header first token)"
    )]
    pub only_id: bool,
    #[arg(
        short = 'b',
        long = "qual-ascii-base",
        default_value_t = 33,
        help = "ASCII base for qualities (usually 33)"
    )]
    pub qual_ascii_base: u8,
    #[arg(
        short = 'g',
        long = "remove-gaps",
        action = ArgAction::SetTrue,
        help = "Remove gap letters configured by -G/--gap-letters"
    )]
    pub remove_gaps: bool,
    #[arg(
        short = 's',
        long = "seq",
        action = ArgAction::SetTrue,
        help = "Only print sequences"
    )]
    pub seq: bool,
    #[arg(
        short = 'v',
        long = "validate-seq",
        action = ArgAction::SetTrue,
        help = "Validate sequence characters against IUPAC nucleotide codes"
    )]
    pub validate_seq: bool,
}

#[derive(Debug, clap::Args)]
#[command(
    after_help = "Examples:\n  seqit fq2fa reads.fq -o reads.fa\n  seqit fq2fa reads.fq.gz -o reads.fa.gz --compression gz"
)]
pub struct Fq2faArgs {
    #[command(flatten)]
    pub io: CommonIoArgs,
}

#[derive(Debug, Clone, ValueEnum)]
pub enum SearchBy {
    Id,
    Name,
    Seq,
    Qual,
}

#[derive(Debug, clap::Args)]
#[command(
    after_help = "Examples:\n  seqit grep reads.fq -p ACTG --by seq -o matched.fq\n  seqit grep reads.fq --pattern-file ids.txt --by id --invert\n  seqit grep --in1 r1.fq --in2 r2.fq -p adapter --pair-mode both -o out.r1.fq -O out.r2.fq\n\nMode notes:\n  - Single-end: provide positional INPUT (or stdin).\n  - Paired-end: provide both --in1 and --in2 together."
)]
pub struct GrepArgs {
    #[command(flatten)]
    pub io: CommonIoArgs,
    #[arg(
        short = '1',
        long = "in1",
        help = "Read 1 input file for paired-end mode"
    )]
    pub in1: Option<String>,
    #[arg(
        short = '2',
        long = "in2",
        help = "Read 2 input file for paired-end mode"
    )]
    pub in2: Option<String>,
    #[arg(
        short = 'O',
        long = "output2",
        help = "Read 2 output file for paired-end mode"
    )]
    pub output2: Option<String>,
    #[arg(short = 'p', long = "pattern", help = "Pattern to match")]
    pub pattern: Option<String>,
    #[arg(
        short = 'f',
        long = "pattern-file",
        help = "File containing one pattern per line"
    )]
    pub pattern_file: Option<String>,
    #[arg(
        short = 'b',
        long = "by",
        value_enum,
        default_value = "id",
        help = "Field to search"
    )]
    pub by: SearchBy,
    #[arg(short = 'x', long = "regex", action = ArgAction::SetTrue, help = "Interpret pattern as regex")]
    pub regex: bool,
    #[arg(short = 'i', long = "ignore-case", action = ArgAction::SetTrue, help = "Case-insensitive matching")]
    pub ignore_case: bool,
    #[arg(short = 'v', long = "invert", action = ArgAction::SetTrue, help = "Select non-matching records")]
    pub invert: bool,
    #[arg(short = 'c', long = "count", action = ArgAction::SetTrue, help = "Print count only")]
    pub count: bool,
    #[arg(short = 'n', long = "only-names", action = ArgAction::SetTrue, help = "Print matching record names only")]
    pub only_names: bool,
    #[arg(short = 'P', long = "pair-mode", default_value = "any", value_parser = ["any", "both"], help = "Pair selection mode")]
    pub pair_mode: String,
    #[arg(
        long = "progress",
        action = ArgAction::SetTrue,
        help = "Show progress on stderr while scanning records"
    )]
    pub progress: bool,
    #[arg(
        long = "allow-unpaired",
        action = ArgAction::SetTrue,
        help = "Continue when paired reads are invalid; skip unpaired records"
    )]
    pub allow_unpaired: bool,
}

#[derive(Debug, clap::Args)]
#[command(
    after_help = "Examples:\n  seqit locate reads.fa -p ACGT\n  seqit locate reads.fa -p 'A[CT]G' --regex --all --bed\n  seqit locate reads.fa --pattern-file motifs.txt --ignore-case"
)]
pub struct LocateArgs {
    #[command(flatten)]
    pub io: CommonIoArgs,
    #[arg(short = 'p', long = "pattern", help = "Pattern to locate")]
    pub pattern: Option<String>,
    #[arg(
        short = 'f',
        long = "pattern-file",
        help = "File containing one pattern per line"
    )]
    pub pattern_file: Option<String>,
    #[arg(short = 'x', long = "regex", action = ArgAction::SetTrue, help = "Interpret pattern as regex")]
    pub regex: bool,
    #[arg(short = 'i', long = "ignore-case", action = ArgAction::SetTrue, help = "Case-insensitive matching")]
    pub ignore_case: bool,
    #[arg(short = 'a', long = "all", action = ArgAction::SetTrue, help = "Report all matches per record")]
    pub all: bool,
    #[arg(short = 'B', long = "bed", action = ArgAction::SetTrue, help = "Emit BED-like coordinates")]
    pub bed: bool,
}

#[derive(Debug, clap::Args)]
#[command(
    after_help = "Examples:\n  seqit sample reads.fq -n 100000 -s 42 -o sub.fq\n  seqit sample reads.fq -r 0.1 -s 42 -o sub.fq\n  seqit sample -i r1.fq -I r2.fq -n 50000 -o sub.r1.fq -O sub.r2.fq\n\nMode notes:\n  - Single-end: provide positional INPUT (or stdin).\n  - Paired-end: provide both mate flags together."
)]
pub struct SampleArgs {
    #[command(flatten)]
    pub io: CommonIoArgs,
    #[arg(
        short = 'i',
        visible_short_alias = '1',
        long = "input",
        help = "Read 1 input file for paired-end mode"
    )]
    pub input1: Option<String>,
    #[arg(
        short = 'I',
        visible_short_alias = '2',
        long = "input2",
        help = "Read 2 input file for paired-end mode"
    )]
    pub input2: Option<String>,
    #[arg(long = "in1", help = "Read 1 input file for paired-end mode")]
    pub in1: Option<String>,
    #[arg(long = "in2", help = "Read 2 input file for paired-end mode")]
    pub in2: Option<String>,
    #[arg(
        short = 'O',
        long = "output2",
        help = "Read 2 output file for paired-end mode"
    )]
    pub output2: Option<String>,
    #[arg(short = 'n', long = "num", help = "Number of records/pairs to sample")]
    pub num: Option<usize>,
    #[arg(short = 'r', long = "rate", help = "Sampling rate in [0,1]")]
    pub rate: Option<f64>,
    #[arg(short = 's', long = "seed", default_value_t = 42, help = "Random seed")]
    pub seed: u64,
    #[arg(
        long = "progress",
        action = ArgAction::SetTrue,
        help = "Show progress on stderr while scanning records"
    )]
    pub progress: bool,
    #[arg(
        long = "allow-unpaired",
        action = ArgAction::SetTrue,
        help = "Continue when paired reads are invalid; skip unpaired records"
    )]
    pub allow_unpaired: bool,
}

#[derive(Debug, Clone, ValueEnum)]
pub enum DupBy {
    Id,
    Seq,
    Full,
}

#[derive(Debug, clap::Args)]
#[command(
    after_help = "Examples:\n  seqit rmdup reads.fa --by seq --keep-first -o dedup.fa\n  seqit rmdup reads.fa --by full --keep-last -o dedup.fa\n  seqit rmdup -i r1.fq -I r2.fq --mark-dup -o d1.fq -O d2.fq\n\nMode notes:\n  - Single-end: provide positional INPUT (or stdin).\n  - Paired-end: provide both mate flags together."
)]
pub struct RmdupArgs {
    #[command(flatten)]
    pub io: CommonIoArgs,
    #[arg(
        short = 'i',
        visible_short_alias = '1',
        long = "input",
        help = "Read 1 input file for paired-end mode"
    )]
    pub input1: Option<String>,
    #[arg(
        short = 'I',
        visible_short_alias = '2',
        long = "input2",
        help = "Read 2 input file for paired-end mode"
    )]
    pub input2: Option<String>,
    #[arg(long = "in1", help = "Read 1 input file for paired-end mode")]
    pub in1: Option<String>,
    #[arg(long = "in2", help = "Read 2 input file for paired-end mode")]
    pub in2: Option<String>,
    #[arg(
        short = 'O',
        long = "output2",
        help = "Read 2 output file for paired-end mode"
    )]
    pub output2: Option<String>,
    #[arg(
        short = 'b',
        long = "by",
        value_enum,
        default_value = "full",
        help = "Duplicate key"
    )]
    pub by: DupBy,
    #[arg(short = 'f', long = "keep-first", action = ArgAction::SetTrue, help = "Keep first occurrence of duplicate")]
    pub keep_first: bool,
    #[arg(short = 'l', long = "keep-last", action = ArgAction::SetTrue, help = "Keep last occurrence of duplicate")]
    pub keep_last: bool,
    #[arg(short = 'c', long = "count-dup", action = ArgAction::SetTrue, help = "Append duplicate count to ID")]
    pub count_dup: bool,
    #[arg(short = 'm', long = "mark-dup", action = ArgAction::SetTrue, help = "Mark duplicate records in ID")]
    pub mark_dup: bool,
    #[arg(
        long = "allow-unpaired",
        action = ArgAction::SetTrue,
        help = "Continue when paired reads are invalid; skip unpaired records"
    )]
    pub allow_unpaired: bool,
}

#[derive(Debug, clap::Args)]
#[command(
    after_help = "Examples:\n  seqit rename reads.fa --prefix sample_ --start 1 --width 8 -o renamed.fa\n      # Generate IDs: sample_000001, sample_000002, ...\n\n  seqit rename reads.fa -e 'lib_PREFIX_N' -p r -w 4 -o renamed.fa\n      # Template supports prefix and index placeholders (see README for exact syntax)\n\n  seqit rename reads.fa --map-file id_map.tsv -o renamed.fa\n      # Exact ID mapping from a 2-column TSV: <old_id>\\t<new_id>\n\n  seqit rename reads.fa --match-regex '^sample_(\\\\d+)$' --replace 'S$1' -o renamed.fa\n      # Regex-based ID rewrite with capture groups\n\n  seqit rename -i r1.fq -I r2.fq --prefix pair_ --keep-pair-suffix -o r1.new.fq -O r2.new.fq\n      # Paired-end renaming while preserving /1 and /2 suffixes\n\nMode notes:\n  - Single-end: provide positional INPUT (or stdin).\n  - Paired-end: provide both mate flags together."
)]
pub struct RenameArgs {
    #[command(flatten)]
    pub io: CommonIoArgs,
    #[arg(
        short = 'i',
        visible_short_alias = '1',
        long = "input",
        help = "Read 1 input file for paired-end mode"
    )]
    pub input1: Option<String>,
    #[arg(
        short = 'I',
        visible_short_alias = '2',
        long = "input2",
        help = "Read 2 input file for paired-end mode"
    )]
    pub input2: Option<String>,
    #[arg(long = "in1", help = "Read 1 input file for paired-end mode")]
    pub in1: Option<String>,
    #[arg(long = "in2", help = "Read 2 input file for paired-end mode")]
    pub in2: Option<String>,
    #[arg(
        short = 'O',
        long = "output2",
        help = "Read 2 output file for paired-end mode"
    )]
    pub output2: Option<String>,
    #[arg(
        short = 'p',
        long = "prefix",
        default_value = "seq",
        help = "Prefix for generated IDs"
    )]
    pub prefix: String,
    #[arg(
        short = 's',
        long = "start",
        default_value_t = 1,
        help = "Starting index"
    )]
    pub start: usize,
    #[arg(
        short = 'w',
        long = "width",
        default_value_t = 6,
        help = "Zero-pad width for index"
    )]
    pub width: usize,
    #[arg(
        short = 'e',
        long = "template",
        help = "Custom template for generated IDs (supports prefix/index placeholders; see README)"
    )]
    pub template: Option<String>,
    #[arg(
        short = 'm',
        long = "map-file",
        value_name = "TSV",
        help = "Tab-separated mapping file: <old_id>\\t<new_id> (no header)"
    )]
    pub map_file: Option<String>,
    #[arg(
        short = 'x',
        long = "match-regex",
        value_name = "PATTERN",
        help = "Regex pattern to match IDs for replacement mode"
    )]
    pub match_regex: Option<String>,
    #[arg(
        short = 'r',
        long = "replace",
        value_name = "REPLACEMENT",
        help = "Replacement text for --match-regex (supports $1, $2, ... groups)"
    )]
    pub replace: Option<String>,
    #[arg(
        short = 'M',
        long = "mode",
        value_enum,
        default_value = "auto",
        help = "Rename mode: auto-detect, sequential generation, mapping, or regex replacement"
    )]
    pub mode: RenameMode,
    #[arg(short = 'k', long = "keep-pair-suffix", action = ArgAction::SetTrue, help = "Preserve /1 and /2 suffixes")]
    pub keep_pair_suffix: bool,
    #[arg(
        long = "allow-unpaired",
        action = ArgAction::SetTrue,
        help = "Continue when paired reads are invalid; skip unpaired records"
    )]
    pub allow_unpaired: bool,
}

#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum RenameMode {
    Auto,
    Generate,
    Map,
    Regex,
}

#[derive(Debug, Clone, ValueEnum)]
pub enum SortBy {
    Id,
    Name,
    Len,
    Seq,
}

#[derive(Debug, clap::Args)]
#[command(
    after_help = "Examples:\n  seqit sort reads.fa --by id -o sorted.fa\n  seqit sort reads.fa --by len --numeric --reverse -o sorted.fa\n  seqit sort reads.fa --by seq --mem 512M --tmp-dir /tmp"
)]
pub struct SortArgs {
    #[command(flatten)]
    pub io: CommonIoArgs,
    #[arg(
        short = 'b',
        long = "by",
        value_enum,
        default_value = "id",
        help = "Sort key"
    )]
    pub by: SortBy,
    #[arg(short = 'n', long = "numeric", action = ArgAction::SetTrue, help = "Use numeric sort where applicable")]
    pub numeric: bool,
    #[arg(short = 'r', long = "reverse", action = ArgAction::SetTrue, help = "Reverse sort order")]
    pub reverse: bool,
    #[arg(
        short = 'T',
        long = "tmp-dir",
        help = "Temporary directory for external sort"
    )]
    pub tmp_dir: Option<String>,
    #[arg(
        short = 'm',
        long = "mem",
        default_value = "128M",
        help = "Memory budget before spilling to disk"
    )]
    pub mem: String,
}

#[derive(Debug, clap::Args)]
#[command(
    after_help = "Examples:\n  seqit shuffle reads.fa -s 42 -o shuffled.fa\n  seqit shuffle reads.fq --mem 1G --tmp-dir /tmp -o shuffled.fq\n  seqit shuffle -i r1.fq -I r2.fq -s 42 -o shuf.r1.fq -O shuf.r2.fq\n\nMode notes:\n  - Single-end: provide positional INPUT (or stdin).\n  - Paired-end: provide both mate flags together."
)]
pub struct ShuffleArgs {
    #[command(flatten)]
    pub io: CommonIoArgs,
    #[arg(
        short = 'i',
        visible_short_alias = '1',
        long = "input",
        help = "Read 1 input file for paired-end mode"
    )]
    pub input1: Option<String>,
    #[arg(
        short = 'I',
        visible_short_alias = '2',
        long = "input2",
        help = "Read 2 input file for paired-end mode"
    )]
    pub input2: Option<String>,
    #[arg(long = "in1", help = "Read 1 input file for paired-end mode")]
    pub in1: Option<String>,
    #[arg(long = "in2", help = "Read 2 input file for paired-end mode")]
    pub in2: Option<String>,
    #[arg(
        short = 'O',
        long = "output2",
        help = "Read 2 output file for paired-end mode"
    )]
    pub output2: Option<String>,
    #[arg(short = 's', long = "seed", default_value_t = 42, help = "Random seed")]
    pub seed: u64,
    #[arg(
        short = 'T',
        long = "tmp-dir",
        help = "Temporary directory for shuffle spill files"
    )]
    pub tmp_dir: Option<String>,
    #[arg(
        short = 'm',
        long = "mem",
        default_value = "128M",
        help = "Memory budget before spilling to disk"
    )]
    pub mem: String,
    #[arg(
        long = "allow-unpaired",
        action = ArgAction::SetTrue,
        help = "Continue when paired reads are invalid; skip unpaired records"
    )]
    pub allow_unpaired: bool,
}

#[derive(Debug, clap::Args)]
#[command(
    after_help = "Examples:\n  seqit spike target.fa -a inserts.fa -s 123 -o spiked.fa\n  seqit spike target.fq -a add.fq -o spiked.fq\n  seqit spike -i t.r1.fq -I t.r2.fq -a a.r1.fq -A a.r2.fq -o out.r1.fq -O out.r2.fq\n\nMode notes:\n  - Single-end: positional INPUT (or stdin), with --add.\n  - Paired-end: provide both input mates and both add mates."
)]
pub struct SpikeArgs {
    #[arg(
        value_name = "INPUT",
        help = "Main input file for single-end mode (default: stdin)"
    )]
    pub input: Option<String>,
    #[arg(
        long = "input",
        hide = true,
        help = "Deprecated alias for single-end main input; prefer positional INPUT"
    )]
    pub input_legacy: Option<String>,
    #[arg(
        short = 'i',
        visible_short_alias = '1',
        long = "input1",
        help = "Read 1 input file for paired-end mode"
    )]
    pub input1: Option<String>,
    #[arg(
        short = 'I',
        long = "input2",
        help = "Read 2 input file for paired-end mode"
    )]
    pub input2: Option<String>,
    #[arg(short = 'a', long = "add", help = "Spike-in file (required)")]
    pub add: String,
    #[arg(
        short = 'A',
        long = "add2",
        help = "Read 2 spike-in file for paired-end mode"
    )]
    pub add2: Option<String>,
    #[arg(long = "in1", help = "Read 1 input file for paired-end mode")]
    pub in1: Option<String>,
    #[arg(long = "in2", help = "Read 2 input file for paired-end mode")]
    pub in2: Option<String>,
    #[arg(
        short = 'o',
        long = "output",
        default_value = "-",
        help = "Output file ('-' for stdout)"
    )]
    pub output: String,
    #[arg(
        short = 'O',
        long = "output2",
        help = "Read 2 output file for paired-end mode"
    )]
    pub output2: Option<String>,
    #[arg(short = 's', long = "seed", default_value_t = 42, help = "Random seed")]
    pub seed: u64,
    #[arg(
        short = 'f',
        long = "format",
        value_enum,
        default_value = "auto",
        help = "Input/output format"
    )]
    pub format: FormatArg,
    #[arg(
        short = 'z',
        long = "compression",
        value_enum,
        default_value = "auto",
        help = "Output compression codec"
    )]
    pub compression: CompressionArg,
    #[arg(
        short = 't',
        long = "threads",
        help = "Number of worker threads to use"
    )]
    pub threads: Option<usize>,
    #[arg(
        long = "allow-unpaired",
        action = ArgAction::SetTrue,
        help = "Continue when paired reads are invalid; skip unpaired records"
    )]
    pub allow_unpaired: bool,
}

#[derive(Debug, clap::Args)]
#[command(
    after_help = "Examples:\n  seqit head reads.fa -n 100 -o first100.fa\n  seqit head reads.fq -p 0.1 -o first10pct.fq\n  seqit head -i r1.fq -I r2.fq -n 50000 -o h1.fq -O h2.fq\n\nMode notes:\n  - Single-end: provide positional INPUT (or stdin).\n  - Paired-end: provide both mate flags together."
)]
pub struct HeadArgs {
    #[command(flatten)]
    pub io: CommonIoArgs,
    #[arg(
        short = 'i',
        visible_short_alias = '1',
        long = "input",
        help = "Read 1 input file for paired-end mode"
    )]
    pub input1: Option<String>,
    #[arg(
        short = 'I',
        visible_short_alias = '2',
        long = "input2",
        help = "Read 2 input file for paired-end mode"
    )]
    pub input2: Option<String>,
    #[arg(
        short = 'O',
        long = "output2",
        help = "Read 2 output file for paired-end mode"
    )]
    pub output2: Option<String>,
    #[arg(short = 'n', long = "num", help = "Number of records/pairs to keep")]
    pub num: Option<usize>,
    #[arg(
        short = 'p',
        long = "proportion",
        help = "Proportion of records/pairs to keep in [0,1]"
    )]
    pub proportion: Option<f64>,
    #[arg(
        long = "allow-unpaired",
        action = ArgAction::SetTrue,
        help = "Continue when paired reads are invalid; skip unpaired records"
    )]
    pub allow_unpaired: bool,
}

#[derive(Debug, clap::Args)]
#[command(
    after_help = "Examples:\n  seqit tail reads.fa -n 100 -o last100.fa\n  seqit tail reads.fq -p 0.1 -o last10pct.fq\n  seqit tail -i r1.fq -I r2.fq -n 50000 -o t1.fq -O t2.fq\n\nMode notes:\n  - Single-end: provide positional INPUT (or stdin).\n  - Paired-end: provide both mate flags together."
)]
pub struct TailArgs {
    #[command(flatten)]
    pub io: CommonIoArgs,
    #[arg(
        short = 'i',
        visible_short_alias = '1',
        long = "input",
        help = "Read 1 input file for paired-end mode"
    )]
    pub input1: Option<String>,
    #[arg(
        short = 'I',
        visible_short_alias = '2',
        long = "input2",
        help = "Read 2 input file for paired-end mode"
    )]
    pub input2: Option<String>,
    #[arg(
        short = 'O',
        long = "output2",
        help = "Read 2 output file for paired-end mode"
    )]
    pub output2: Option<String>,
    #[arg(short = 'n', long = "num", help = "Number of records/pairs to keep")]
    pub num: Option<usize>,
    #[arg(
        short = 'p',
        long = "proportion",
        help = "Proportion of records/pairs to keep in [0,1]"
    )]
    pub proportion: Option<f64>,
    #[arg(
        long = "allow-unpaired",
        action = ArgAction::SetTrue,
        help = "Continue when paired reads are invalid; skip unpaired records"
    )]
    pub allow_unpaired: bool,
}
