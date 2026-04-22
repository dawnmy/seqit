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
        }
    }
}

#[derive(Debug, clap::Args, Clone)]
pub struct CommonIoArgs {
    #[arg(value_name = "INPUT")]
    pub input: Option<String>,
    #[arg(short = 'o', long = "output", default_value = "-")]
    pub output: String,
    #[arg(long = "format", value_enum, default_value = "auto")]
    pub format: FormatArg,
    #[arg(long = "compression", value_enum, default_value = "auto")]
    pub compression: CompressionArg,
    #[arg(
        short = 't',
        long = "threads",
        help = "Number of worker threads to use"
    )]
    pub threads: Option<usize>,
    #[arg(long = "quiet", action = ArgAction::SetTrue)]
    pub quiet: bool,
}

#[derive(Debug, clap::Args)]
pub struct StatsArgs {
    #[arg(value_name = "INPUTS")]
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
        short = 'p',
        long = "per-file",
        action = ArgAction::SetTrue,
        help = "When multiple files are given, do not append a TOTAL summary row"
    )]
    pub per_file: bool,
    #[arg(long = "format", value_enum, default_value = "auto")]
    pub format: FormatArg,
    #[arg(
        short = 't',
        long = "threads",
        help = "Number of worker threads to use"
    )]
    pub threads: Option<usize>,
}

#[derive(Debug, clap::Args)]
pub struct SeqArgs {
    #[command(flatten)]
    pub io: CommonIoArgs,
    #[arg(short = 'm', long = "min-len")]
    pub min_len: Option<usize>,
    #[arg(short = 'M', long = "max-len")]
    pub max_len: Option<usize>,
    #[arg(short = 'r', long = "rev", action = ArgAction::SetTrue)]
    pub rev: bool,
    #[arg(short = 'c', long = "comp", action = ArgAction::SetTrue)]
    pub comp: bool,
    #[arg(short = 'R', long = "revcomp", action = ArgAction::SetTrue)]
    pub revcomp: bool,
    #[arg(short = 'u', long = "upper", action = ArgAction::SetTrue)]
    pub upper: bool,
    #[arg(short = 'l', long = "lower", action = ArgAction::SetTrue)]
    pub lower: bool,
}

#[derive(Debug, clap::Args)]
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
pub struct GrepArgs {
    #[command(flatten)]
    pub io: CommonIoArgs,
    #[arg(short = '1', long = "in1")]
    pub in1: Option<String>,
    #[arg(short = '2', long = "in2")]
    pub in2: Option<String>,
    #[arg(short = 'O', long = "output2")]
    pub output2: Option<String>,
    #[arg(short = 'p', long = "pattern")]
    pub pattern: Option<String>,
    #[arg(short = 'f', long = "pattern-file")]
    pub pattern_file: Option<String>,
    #[arg(short = 'b', long = "by", value_enum, default_value = "id")]
    pub by: SearchBy,
    #[arg(short = 'x', long = "regex", action = ArgAction::SetTrue)]
    pub regex: bool,
    #[arg(short = 'i', long = "ignore-case", action = ArgAction::SetTrue)]
    pub ignore_case: bool,
    #[arg(short = 'v', long = "invert", action = ArgAction::SetTrue)]
    pub invert: bool,
    #[arg(short = 'c', long = "count", action = ArgAction::SetTrue)]
    pub count: bool,
    #[arg(short = 'n', long = "only-names", action = ArgAction::SetTrue)]
    pub only_names: bool,
    #[arg(short = 'P', long = "pair-mode", default_value = "any", value_parser = ["any", "both"])]
    pub pair_mode: String,
}

#[derive(Debug, clap::Args)]
pub struct LocateArgs {
    #[command(flatten)]
    pub io: CommonIoArgs,
    #[arg(short = 'p', long = "pattern")]
    pub pattern: Option<String>,
    #[arg(short = 'f', long = "pattern-file")]
    pub pattern_file: Option<String>,
    #[arg(short = 'x', long = "regex", action = ArgAction::SetTrue)]
    pub regex: bool,
    #[arg(short = 'i', long = "ignore-case", action = ArgAction::SetTrue)]
    pub ignore_case: bool,
    #[arg(short = 'a', long = "all", action = ArgAction::SetTrue)]
    pub all: bool,
    #[arg(short = 'B', long = "bed", action = ArgAction::SetTrue)]
    pub bed: bool,
}

#[derive(Debug, clap::Args)]
pub struct SampleArgs {
    #[command(flatten)]
    pub io: CommonIoArgs,
    #[arg(short = '1', long = "in1")]
    pub in1: Option<String>,
    #[arg(short = '2', long = "in2")]
    pub in2: Option<String>,
    #[arg(short = 'O', long = "output2")]
    pub output2: Option<String>,
    #[arg(short = 'n', long = "num")]
    pub num: Option<usize>,
    #[arg(short = 'r', long = "rate")]
    pub rate: Option<f64>,
    #[arg(short = 's', long = "seed", default_value_t = 42)]
    pub seed: u64,
}

#[derive(Debug, Clone, ValueEnum)]
pub enum DupBy {
    Id,
    Seq,
    Full,
}

#[derive(Debug, clap::Args)]
pub struct RmdupArgs {
    #[command(flatten)]
    pub io: CommonIoArgs,
    #[arg(short = '1', long = "in1")]
    pub in1: Option<String>,
    #[arg(short = '2', long = "in2")]
    pub in2: Option<String>,
    #[arg(short = 'O', long = "output2")]
    pub output2: Option<String>,
    #[arg(short = 'b', long = "by", value_enum, default_value = "full")]
    pub by: DupBy,
    #[arg(short = 'f', long = "keep-first", action = ArgAction::SetTrue)]
    pub keep_first: bool,
    #[arg(short = 'l', long = "keep-last", action = ArgAction::SetTrue)]
    pub keep_last: bool,
    #[arg(short = 'c', long = "count-dup", action = ArgAction::SetTrue)]
    pub count_dup: bool,
    #[arg(short = 'm', long = "mark-dup", action = ArgAction::SetTrue)]
    pub mark_dup: bool,
}

#[derive(Debug, clap::Args)]
pub struct RenameArgs {
    #[command(flatten)]
    pub io: CommonIoArgs,
    #[arg(short = '1', long = "in1")]
    pub in1: Option<String>,
    #[arg(short = '2', long = "in2")]
    pub in2: Option<String>,
    #[arg(short = 'O', long = "output2")]
    pub output2: Option<String>,
    #[arg(short = 'p', long = "prefix", default_value = "seq")]
    pub prefix: String,
    #[arg(short = 's', long = "start", default_value_t = 1)]
    pub start: usize,
    #[arg(short = 'w', long = "width", default_value_t = 6)]
    pub width: usize,
    #[arg(short = 'e', long = "template")]
    pub template: Option<String>,
    #[arg(short = 'k', long = "keep-pair-suffix", action = ArgAction::SetTrue)]
    pub keep_pair_suffix: bool,
}

#[derive(Debug, Clone, ValueEnum)]
pub enum SortBy {
    Id,
    Name,
    Len,
    Seq,
}

#[derive(Debug, clap::Args)]
pub struct SortArgs {
    #[command(flatten)]
    pub io: CommonIoArgs,
    #[arg(short = 'b', long = "by", value_enum, default_value = "id")]
    pub by: SortBy,
    #[arg(short = 'n', long = "numeric", action = ArgAction::SetTrue)]
    pub numeric: bool,
    #[arg(short = 'r', long = "reverse", action = ArgAction::SetTrue)]
    pub reverse: bool,
    #[arg(short = 'T', long = "tmp-dir")]
    pub tmp_dir: Option<String>,
    #[arg(short = 'm', long = "mem", default_value = "128M")]
    pub mem: String,
}

#[derive(Debug, clap::Args)]
pub struct ShuffleArgs {
    #[command(flatten)]
    pub io: CommonIoArgs,
    #[arg(short = '1', long = "in1")]
    pub in1: Option<String>,
    #[arg(short = '2', long = "in2")]
    pub in2: Option<String>,
    #[arg(short = 'O', long = "output2")]
    pub output2: Option<String>,
    #[arg(short = 's', long = "seed", default_value_t = 42)]
    pub seed: u64,
    #[arg(short = 'T', long = "tmp-dir")]
    pub tmp_dir: Option<String>,
    #[arg(short = 'm', long = "mem", default_value = "128M")]
    pub mem: String,
}

#[derive(Debug, clap::Args)]
pub struct SpikeArgs {
    #[arg(short = 'i', long = "input")]
    pub input: Option<String>,
    #[arg(short = 'a', long = "add")]
    pub add: String,
    #[arg(short = 'A', long = "add2")]
    pub add2: Option<String>,
    #[arg(short = '1', long = "in1")]
    pub in1: Option<String>,
    #[arg(short = '2', long = "in2")]
    pub in2: Option<String>,
    #[arg(short = 'o', long = "output", default_value = "-")]
    pub output: String,
    #[arg(short = 'O', long = "output2")]
    pub output2: Option<String>,
    #[arg(short = 's', long = "seed", default_value_t = 42)]
    pub seed: u64,
    #[arg(short = 'f', long = "format", value_enum, default_value = "auto")]
    pub format: FormatArg,
    #[arg(short = 'z', long = "compression", value_enum, default_value = "auto")]
    pub compression: CompressionArg,
    #[arg(
        short = 't',
        long = "threads",
        help = "Number of worker threads to use"
    )]
    pub threads: Option<usize>,
}
