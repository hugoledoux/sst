//-- sstsimpl

#[allow(dead_code)]
#[allow(unused_variables)]
#[macro_use]
extern crate log; //info/debug/error

extern crate startin;

use std::io::BufRead;
use std::io::{self, Write};

// use std::cmp::Ordering;
// use std::collections::BinaryHeap;

use clap::Parser;
#[derive(Parser)]
#[command(name = "sstsimpl")]
#[command(about = "streaming startin -- simplify the terrain [sstsimpl]")]
#[command(author, version)]
struct Cli {
    /// vertical epsilon
    vepsilon: f64,
    #[clap(flatten)]
    verbose: clap_verbosity_flag::Verbosity,
}

// struct PotentialVertex {
//     id: usize,
//     verr: f64,
// }

// impl PartialOrd for PotentialVertex {
//     fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
//         self.verr.partial_cmp(&other.verr)
//     }
// }

// impl PartialEq for PotentialVertex {
//     fn eq(&self, other: &Self) -> bool {
//         self.verr == other.verr
//     }
// }

struct Surface {
    pts: Vec<Vec<f64>>,
    max_epsilon: f64,
}

impl Surface {
    fn new(epsilon: f64) -> Surface {
        let l = Vec::new();
        Surface {
            pts: l,
            max_epsilon: epsilon,
        }
    }

    fn clear(&mut self) {
        self.pts.clear();
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
            if p[1] > ymax {
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

    fn finalise(&mut self) -> io::Result<()> {
        info!("finalise: {}", self.pts.len());
        info!("max_epsilon: {}", self.max_epsilon);
        let bbox = self.get_bbox();
        let zavg = self.get_average_elevation();
        let mut dt = startin::Triangulation::new();
        //-- insert 4 dummy corner TODO: better z-values would be nice
        let _ = dt.insert_one_pt(bbox[0], bbox[1], zavg);
        let _ = dt.insert_one_pt(bbox[2], bbox[1], zavg);
        let _ = dt.insert_one_pt(bbox[2], bbox[3], zavg);
        let _ = dt.insert_one_pt(bbox[0], bbox[3], zavg);
        let mut epsilon = std::f64::MAX;
        while epsilon > self.max_epsilon {
            // info!("self.pts.len()= {}", self.pts.len());
            epsilon = 0.0;
            let mut pi: usize = 0;
            for (i, p) in self.pts.iter().enumerate() {
                let z2 = dt.interpolate_tin_linear(p[0], p[1]).unwrap();
                let e = (p[2] - z2).abs();
                if e > epsilon {
                    // info!("{:?}", e);
                    epsilon = e;
                    pi = i;
                }
            }
            let _ = dt.insert_one_pt(self.pts[pi][0], self.pts[pi][1], self.pts[pi][2]);
            self.pts.remove(pi);
        }

        let impdigits = 3;
        let allv = &dt.all_vertices();
        for i in 5..dt.number_of_vertices() {
            io::stdout().write_all(
                &format!(
                    "v {0:.3$} {1:.3$} {2:.3$}\n",
                    allv[i][0], allv[i][1], allv[i][2], impdigits
                )
                .as_bytes(),
            )?;
        }
        Ok(())
    }
}

fn main() -> io::Result<()> {
    let cli = Cli::parse();
    env_logger::Builder::new()
        .filter_level(cli.verbose.log_level_filter())
        .init();

    // let mut pts: Vec<Vec<f64>> = Vec::new();
    let mut s = Surface::new(cli.vepsilon);

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
                //-- finalise a cell
                if s.is_empty() {
                    io::stdout().write_all(&format!("{}\n", &l).as_bytes())?;
                } else {
                    let _ = s.finalise();
                    s.clear();
                    io::stdout().write_all(&format!("{}\n", &l).as_bytes())?;
                }
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

fn parse_3_f64(l: &String) -> Vec<f64> {
    let ls: Vec<&str> = l.split_whitespace().collect();
    let a: f64 = ls[1].parse::<f64>().unwrap();
    let b: f64 = ls[2].parse::<f64>().unwrap();
    let c: f64 = ls[3].parse::<f64>().unwrap();
    vec![a, b, c]
}
