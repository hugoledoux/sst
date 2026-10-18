//-- ttdt

#[allow(dead_code)]
#[allow(unused_variables)]
mod triangulator;

mod gui;

#[macro_use]
extern crate log; //info/debug/error

use std::io::BufRead;
use std::io::{self, Write};

use clap::Parser;

use triangulator::Triangulation;

#[derive(Parser)]
#[command(name = "ttdt")]
#[command(about = "tintin -- Delaunay triangulator [ttdt]")]
#[command(author, version)]
pub(crate) struct Cli {
    /// Show a live 2D view of the current triangulation and the quadtree cells
    #[arg(long)]
    gui: bool,
    #[clap(flatten)]
    verbose: clap_verbosity_flag::Verbosity,
}

fn main() -> io::Result<()> {
    let cli = Cli::parse();
    env_logger::Builder::new()
        .filter_level(cli.verbose.log_level_filter())
        .init();

    if cli.gui {
        gui::run()
    } else {
        run_cli()
    }
}

fn run_cli() -> io::Result<()> {
    info!("Init DT");
    let mut dt = Triangulation::new();

    //----- reading from stdin -----//
    let stdin = std::io::stdin();
    for line in stdin.lock().lines() {
        let l = line.unwrap();
        handle_line(&mut dt, &l)?;
    }

    let _re = dt.flush_output();
    log_stats(&dt);
    Ok(())
}

/// Process one record of the `.spa` stream. Shared by the CLI loop and the GUI.
pub(crate) fn handle_line(dt: &mut Triangulation, l: &str) -> io::Result<()> {
    if l.is_empty() {
        return Ok(());
    }
    let ch = l.chars().next().unwrap();
    match ch {
        '#' => {}
        'n' => {
            //-- number of points (not needed)
            let _totalpts = l.split_whitespace().last().unwrap().parse::<usize>().unwrap();
        }
        's' => {
            //-- cellsize
            let s = l.split_whitespace().last().unwrap().parse::<usize>().unwrap();
            dt.set_cellsize(s);
        }
        'c' => {
            //-- dimension grid (always square cXc)
            let c = l.split_whitespace().last().unwrap().parse::<usize>().unwrap();
            dt.set_grid_dimensions(c);
        }
        'b' => {
            //-- bbox
            let re = parse_4_f64(l);
            dt.set_bbox(re.0, re.1, re.2, re.3);
            io::stdout().write_all(
                &format!("b {:.3} {:.3} {:.3} {:.3}\n", re.0, re.1, re.2, re.3).as_bytes(),
            )?;
        }
        'v' | 'w' => {
            //-- vertex
            let v = parse_3_f64(l);
            let _re = dt.insert_one_pt_with_grid(v.0, v.1, v.2);
        }
        'x' => {
            //-- finalise a cell
            let re = parse_2_usize(l);
            let _re = dt.finalise_qt_leaf(re.0, re.1);
        }
        _ => {
            error!("Wrongly formatted stream. Abort.");
            std::process::exit(1);
        }
    }
    Ok(())
}

fn log_stats(dt: &Triangulation) {
    info!("dt.number_of_vertices() = {}", dt.number_of_vertices());

    info!("max # points in DT during process: {}", dt.max_vs);
    info!("max # triangles in DT during process: {}", dt.max_ts);

    info!(
        "# insertions: {}; # flips22: {} (avg {:.2} per insertion)",
        dt.n_inserts,
        dt.n_flips22,
        dt.n_flips22 as f64 / dt.n_inserts as f64
    );
    info!(
        "# walks: {}; # walk steps: {} (avg {:.2} per walk)",
        dt.n_walks,
        dt.n_walksteps,
        dt.n_walksteps as f64 / dt.n_walks as f64
    );

    info!("✅");
}

fn parse_2_usize(l: &str) -> (usize, usize) {
    let ls: Vec<&str> = l.split_whitespace().collect();
    let a: usize = ls[1].parse::<usize>().unwrap();
    let b: usize = ls[2].parse::<usize>().unwrap();
    (a, b)
}

fn parse_3_f64(l: &str) -> (f64, f64, f64) {
    let ls: Vec<&str> = l.split_whitespace().collect();
    let a: f64 = ls[1].parse::<f64>().unwrap();
    let b: f64 = ls[2].parse::<f64>().unwrap();
    let c: f64 = ls[3].parse::<f64>().unwrap();
    (a, b, c)
}

fn parse_4_f64(l: &str) -> (f64, f64, f64, f64) {
    let ls: Vec<&str> = l.split_whitespace().collect();
    let a: f64 = ls[1].parse::<f64>().unwrap();
    let b: f64 = ls[2].parse::<f64>().unwrap();
    let c: f64 = ls[3].parse::<f64>().unwrap();
    let d: f64 = ls[4].parse::<f64>().unwrap();
    (a, b, c, d)
}
