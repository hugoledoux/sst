//-- sstcount

#[allow(dead_code)]
#[allow(unused_variables)]
#[macro_use]
extern crate log; //info/debug/error

use clap::Parser;
use num_format::{Locale, ToFormattedString};
use std::io::BufRead;
use std::io::{self, Write};

#[derive(Parser)]
#[command(name = "sstcount")]
#[command(about = "streaming startin -- count the vertices/triangles [sstcount]")]
#[command(author, version)]
struct Cli {
    #[clap(flatten)]
    verbose: clap_verbosity_flag::Verbosity,
}

fn main() -> io::Result<()> {
    let cli = Cli::parse();
    env_logger::Builder::new()
        .filter_level(cli.verbose.log_level_filter())
        .init();

    let mut total_vs: usize = 0;
    let mut total_ts: usize = 0;

    let stdin = std::io::stdin();
    for line in stdin.lock().lines() {
        let l = line.unwrap();
        if l.is_empty() {
            continue;
        }
        io::stdout().write_all(&format!("{}\n", &l).as_bytes())?;
        let ch = l.chars().next().unwrap();
        match ch {
            'v' => {
                //-- a vertex
                total_vs += 1;
            }
            'f' => {
                //-- a triangle
                total_ts += 1;
            }
            _ => {}
        }
    }
    info!("# vertices: {}", total_vs.to_formatted_string(&Locale::en));
    info!("# triangles: {}", total_ts.to_formatted_string(&Locale::en));
    info!("✅");
    Ok(())
}
