//-- sstsimpl

#[allow(dead_code)]
#[allow(unused_variables)]
#[macro_use]
extern crate log; //info/debug/error

extern crate startin;

use std::io::BufRead;
use std::io::{self, Write};

use clap::Parser;
#[derive(Parser)]
#[command(name = "sstsimpl")]
#[command(about = "streaming startin -- simplify the terrain [sstsimpl]")]
#[command(author, version)]
struct Cli {
    #[clap(flatten)]
    verbose: clap_verbosity_flag::Verbosity,
}

struct Surface {
    pts: Vec<Vec<f64>>,
}

impl Surface {
    fn new() -> Surface {
        let l = Vec::new();
        Surface { pts: l }
    }

    fn add_pt(&mut self, pt: Vec<f64>) {
        self.pts.push(pt);
    }

    fn is_empty(&self) -> bool {
        self.pts.is_empty()
    }

    fn get_bbox(&self) -> Vec<f64> {
        let mut xmin: f64 = f64::MAX;
        let mut ymin: f64 = f64::MAX;
        let mut xmax: f64 = f64::MIN;
        let mut ymax: f64 = f64::MIN;
        for p in &self.pts {
            if p[0] < xmin {
                xmin = p[0]
            }
            if p[1] < ymin {
                ymin = p[1]
            }
            if p[0] > xmax {
                xmax = p[0]
            }
            if p[1] < ymax {
                ymax = p[1]
            }
        }
        vec![xmin - 10.0, ymin - 10.0, xmax + 10.0, ymax + 10.0]
    }

    fn get_average_elevation(&self) -> f64 {
        let mut total = 0_f64;
        let mut n: usize = 0;
        for p in self.pts.iter().step_by(10) {
            total += p[2] as f64;
            n += 1;
        }
        total / n as f64
    }

    fn finalise(&self) -> io::Result<()> {
        let bbox = self.get_bbox();
        let zavg = self.get_average_elevation();
        let mut dt = startin::Triangulation::new();
        //-- insert 4 dummy corner TODO: better z-values would be nice
        let _ = dt.insert_one_pt(bbox[0], bbox[1], zavg);
        let _ = dt.insert_one_pt(bbox[2], bbox[1], zavg);
        let _ = dt.insert_one_pt(bbox[2], bbox[3], zavg);
        let _ = dt.insert_one_pt(bbox[0], bbox[3], zavg);

        // dt.insert(&self.get_bbox(), startin::InsertionStrategy::AsIs);

        // io::stdout().write_all(&format!("{}\n", "FINALISE").as_bytes())?;
        Ok(())
    }
}

fn main() -> io::Result<()> {
    let cli = Cli::parse();
    env_logger::Builder::new()
        .filter_level(cli.verbose.log_level_filter())
        .init();

    // let mut pts: Vec<Vec<f64>> = Vec::new();
    let mut s = Surface::new();

    let stdin = std::io::stdin();
    for line in stdin.lock().lines() {
        let l = line.unwrap();
        if l.is_empty() {
            continue;
        }

        let ch = l.chars().next().unwrap();
        match ch {
            '#' => continue,
            'n' => {
                //-- number of points
                io::stdout().write_all(&format!("{}\n", &l).as_bytes())?;
            }
            's' => {
                //-- cellsize
                io::stdout().write_all(&format!("{}\n", &l).as_bytes())?;
            }
            'c' => {
                //-- dimension grid (always square cXc)
                io::stdout().write_all(&format!("{}\n", &l).as_bytes())?;
            }
            'b' => {
                //-- bbox
                io::stdout().write_all(&format!("{}\n", &l).as_bytes())?;
            }
            'v' => {
                //-- a vertex
                s.add_pt(parse_3_f64(&l));
            }
            'x' => {
                //-- finalise a vertex
                if s.is_empty() {
                    io::stdout().write_all(&format!("{}\n", &l).as_bytes())?;
                } else {
                    //-- create a DT
                    s.finalise();
                }
                // let ls: Vec<&str> = l.split_whitespace().collect();
                // let a: usize = ls[1].parse::<usize>().unwrap();
                // d.remove(&a);
                continue;
            }
            _ => {
                error!("Wrongly formatted stream. Abort.");
                std::process::exit(1);
            }
        }
    }

    info!("✅");
    Ok(())
}

// fn finalise_cell(pts: &Vec<Vec<f64>>) {
//     let mut dt = startin::Triangulation::new();
//     //-- insert 4 corner
//     let total = 0_f64;
//     for pt in pts {
//         total += pt[2];
//     }
//     let avgz: f64 = total / pts.len();

//     //-- insert all
//     // dt.insert(&pts);
// }

fn parse_3_f64(l: &String) -> Vec<f64> {
    let ls: Vec<&str> = l.split_whitespace().collect();
    let a: f64 = ls[1].parse::<f64>().unwrap();
    let b: f64 = ls[2].parse::<f64>().unwrap();
    let c: f64 = ls[3].parse::<f64>().unwrap();
    vec![a, b, c]
}

fn parse_4_f64(l: &String) -> (f64, f64, f64, f64) {
    let ls: Vec<&str> = l.split_whitespace().collect();
    let a: f64 = ls[1].parse::<f64>().unwrap();
    let b: f64 = ls[2].parse::<f64>().unwrap();
    let c: f64 = ls[3].parse::<f64>().unwrap();
    let d: f64 = ls[4].parse::<f64>().unwrap();
    (a, b, c, d)
}
