//-- spb2obj

#[allow(dead_code)]
#[allow(unused_variables)]
#[macro_use]
extern crate log; //info/debug/error

use std::io::BufRead;
use std::io::{self, Write};

use clap::Parser;
#[derive(Parser)]
#[command(name = "spb2obj")]
#[command(about = "tintin -- convert the .sma stream of ttdt to a Wavefront OBJ [spb2obj]")]
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

    let stdout = io::stdout();
    let mut w = std::io::BufWriter::with_capacity(1 << 20, stdout.lock());

    let stdin = std::io::stdin();
    let mut r = std::io::BufReader::new(stdin.lock());
    if sst::stream::peek_is_binary(&mut r, sst::stream::SMA_MAGIC)? {
        info!("input stream: binary .sma");
        sst::stream::skip_magic(&mut r, sst::stream::SMA_MAGIC)?;
        while let Some(rec) = sst::stream::sma_read(&mut r)? {
            match rec {
                sst::stream::SmaRecord::V(..) | sst::stream::SmaRecord::F(..) => {
                    w.write_all(rec.to_line().as_bytes())?;
                }
                //-- the bbox and the vertex-finalisation records have no
                //-- place in an OBJ file
                sst::stream::SmaRecord::B(..) | sst::stream::SmaRecord::X(..) => {}
            }
        }
    } else if sst::stream::peek_is_binary(&mut r, sst::stream::SPA_MAGIC)? {
        error!("this is a .spa stream (points); spb2obj expects the .sma stream of ttdt");
        std::process::exit(1);
    } else {
        info!("input stream: ASCII .sma");
        for line in r.lines() {
            let l = line.unwrap();
            if l.is_empty() {
                continue;
            }
            let ch = l.chars().next().unwrap();
            match ch {
                'v' | 'f' => {
                    w.write_all(l.as_bytes())?;
                    w.write_all(b"\n")?;
                }
                '#' | 'b' | 'x' => {}
                _ => {
                    error!("Wrongly formatted stream. Abort.");
                    std::process::exit(1);
                }
            }
        }
    }
    w.flush()?;
    info!("✅");
    Ok(())
}
