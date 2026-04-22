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
    Stats(StatsArgs),
    Seq(SeqArgs),
    Fq2fa(Fq2faArgs),
    Grep(GrepArgs),
    Locate(LocateArgs),
    Sample(SampleArgs),
    Rmdup(RmdupArgs),
    Rename(RenameArgs),
    Sort(SortArgs),
    Shuffle(ShuffleArgs),
    Spike(SpikeArgs),
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
    #[arg(short = 't', long = "threads")]
    pub threads: Option<usize>,
    #[arg(long = "quiet", action = ArgAction::SetTrue)]
    pub quiet: bool,
}

#[derive(Debug, clap::Args)]
pub struct StatsArgs {
    #[arg(value_name = "INPUTS")]
    pub inputs: Vec<String>,
    #[arg(long = "tabular", action = ArgAction::SetTrue)]
    pub tabular: bool,
    #[arg(long = "json", action = ArgAction::SetTrue)]
    pub json: bool,
    #[arg(long = "per-file", action = ArgAction::SetTrue)]
    pub per_file: bool,
    #[arg(long = "format", value_enum, default_value = "auto")]
    pub format: FormatArg,
}

#[derive(Debug, clap::Args)]
pub struct SeqArgs {
    #[command(flatten)]
    pub io: CommonIoArgs,
    #[arg(long = "min-len")]
    pub min_len: Option<usize>,
    #[arg(long = "max-len")]
    pub max_len: Option<usize>,
    #[arg(long = "rev", action = ArgAction::SetTrue)]
    pub rev: bool,
    #[arg(long = "comp", action = ArgAction::SetTrue)]
    pub comp: bool,
    #[arg(long = "revcomp", action = ArgAction::SetTrue)]
    pub revcomp: bool,
    #[arg(long = "upper", action = ArgAction::SetTrue)]
    pub upper: bool,
    #[arg(long = "lower", action = ArgAction::SetTrue)]
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
    #[arg(long = "by", value_enum, default_value = "id")]
    pub by: SearchBy,
    #[arg(long = "regex", action = ArgAction::SetTrue)]
    pub regex: bool,
    #[arg(long = "ignore-case", action = ArgAction::SetTrue)]
    pub ignore_case: bool,
    #[arg(short = 'v', long = "invert", action = ArgAction::SetTrue)]
    pub invert: bool,
    #[arg(short = 'c', long = "count", action = ArgAction::SetTrue)]
    pub count: bool,
    #[arg(short = 'n', long = "only-names", action = ArgAction::SetTrue)]
    pub only_names: bool,
    #[arg(long = "pair-mode", default_value = "any", value_parser = ["any", "both"])]
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
    #[arg(long = "regex", action = ArgAction::SetTrue)]
    pub regex: bool,
    #[arg(long = "ignore-case", action = ArgAction::SetTrue)]
    pub ignore_case: bool,
    #[arg(long = "all", action = ArgAction::SetTrue)]
    pub all: bool,
    #[arg(long = "bed", action = ArgAction::SetTrue)]
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
    #[arg(long = "by", value_enum, default_value = "full")]
    pub by: DupBy,
    #[arg(long = "keep-first", action = ArgAction::SetTrue)]
    pub keep_first: bool,
    #[arg(long = "keep-last", action = ArgAction::SetTrue)]
    pub keep_last: bool,
    #[arg(long = "count-dup", action = ArgAction::SetTrue)]
    pub count_dup: bool,
    #[arg(long = "mark-dup", action = ArgAction::SetTrue)]
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
    #[arg(long = "prefix", default_value = "seq")]
    pub prefix: String,
    #[arg(long = "start", default_value_t = 1)]
    pub start: usize,
    #[arg(long = "width", default_value_t = 6)]
    pub width: usize,
    #[arg(long = "template")]
    pub template: Option<String>,
    #[arg(long = "keep-pair-suffix", action = ArgAction::SetTrue)]
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
    #[arg(long = "by", value_enum, default_value = "id")]
    pub by: SortBy,
    #[arg(long = "numeric", action = ArgAction::SetTrue)]
    pub numeric: bool,
    #[arg(long = "reverse", action = ArgAction::SetTrue)]
    pub reverse: bool,
    #[arg(long = "tmp-dir")]
    pub tmp_dir: Option<String>,
    #[arg(long = "mem", default_value = "128M")]
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
    #[arg(long = "tmp-dir")]
    pub tmp_dir: Option<String>,
    #[arg(long = "mem", default_value = "128M")]
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
    #[arg(long = "format", value_enum, default_value = "auto")]
    pub format: FormatArg,
    #[arg(long = "compression", value_enum, default_value = "auto")]
    pub compression: CompressionArg,
}
