mod cli;
mod commands;
mod formats;
mod io;
mod pairs;
mod utils;

use anyhow::Result;
use clap::Parser;
use cli::{Cli, Commands};
use std::io::ErrorKind;

fn main() -> Result<()> {
    let cli = Cli::parse();
    let default_threads = std::thread::available_parallelism()
        .map(|n| n.get().min(10))
        .unwrap_or(1);
    let threads = cli
        .command
        .threads()
        .filter(|&t| t > 0)
        .unwrap_or(default_threads);
    if threads > 0 {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global();
    }
    let res = match cli.command {
        Commands::Stats(a) => commands::stats::run(a),
        Commands::Seq(a) => commands::seq::run(a),
        Commands::Fq2fa(a) => commands::fq2fa::run(a),
        Commands::Grep(a) => commands::grep::run(a),
        Commands::Locate(a) => commands::locate::run(a),
        Commands::Sample(a) => commands::sample::run(a),
        Commands::Rmdup(a) => commands::rmdup::run(a),
        Commands::Rename(a) => commands::rename::run(a),
        Commands::Sort(a) => commands::sort::run(a),
        Commands::Shuffle(a) => commands::shuffle::run(a),
        Commands::Spike(a) => commands::spike::run(a),
        Commands::Head(a) => commands::head::run(a),
        Commands::Tail(a) => commands::tail::run(a),
    };
    if let Err(err) = res {
        if err.chain().any(|cause| {
            cause
                .downcast_ref::<std::io::Error>()
                .is_some_and(|e| e.kind() == ErrorKind::BrokenPipe)
        }) {
            return Ok(());
        }
        return Err(err);
    }
    Ok(())
}
