mod cli;
mod commands;
mod formats;
mod io;
mod pairs;
mod utils;

use anyhow::Result;
use clap::Parser;
use cli::{Cli, Commands};

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
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
    }
}
