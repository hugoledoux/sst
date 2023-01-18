//-- sstsimpl2

#[allow(dead_code)]
#[allow(unused_variables)]
#[macro_use]
extern crate log; //info/debug/error

extern crate startin;

extern crate rayon;
use rayon::prelude::*;

use std::io::BufRead;
use std::io::{self, Write};

use rand::seq::SliceRandom;
use rand::thread_rng;

use clap::Parser;
#[derive(Parser)]
#[command(name = "sstsimpl2")]
#[command(about = "streaming startin -- simplify the terrain [parallel version]")]
#[command(author, version)]
struct Cli {
    /// Vertical epsilon
    vepsilon: f64,
    /// Number of cores to use
    #[arg(short, long, default_value_t = 4)]
    cores: usize,
    #[clap(flatten)]
    verbose: clap_verbosity_flag::Verbosity,
}

#[derive(Clone)]
struct Surface {
    pts: Vec<Vec<f64>>,
    max_epsilon: f64,
    bbox: Vec<f64>,
    cell: String,
}

impl Surface {
    fn new(epsilon: f64) -> Surface {
        let l = Vec::new();
        let b = vec![f64::MAX, f64::MAX, f64::MIN, f64::MIN];
        let s = "x 0 0";
        Surface {
            pts: l,
            max_epsilon: epsilon,
            bbox: b,
            cell: s.to_string(),
        }
    }

    fn add_pt(&mut self, p: Vec<f64>) {
        if p[0] < self.bbox[0] {
            self.bbox[0] = p[0]
        }
        if p[1] < self.bbox[1] {
            self.bbox[1] = p[1]
        }
        if p[0] > self.bbox[2] {
            self.bbox[2] = p[0]
        }
        if p[1] > self.bbox[3] {
            self.bbox[3] = p[1]
        }
        self.pts.push(p);
    }

    fn is_empty(&self) -> bool {
        self.pts.is_empty()
    }

    fn close_to_bbox(&self, p: &Vec<f64>, d: f64) -> bool {
        // let bbox = self.get_bbox();
        let deltax = self.bbox[2] - self.bbox[0];
        let deltay = self.bbox[3] - self.bbox[1];
        if p[0] - self.bbox[0] < (d * deltax) {
            return true;
        }
        if self.bbox[2] - p[0] < (d * deltax) {
            return true;
        }
        if p[1] - self.bbox[1] < (d * deltay) {
            return true;
        }
        if self.bbox[3] - p[1] < (d * deltay) {
            return true;
        }
        false
    }

    fn get_corner_elevations(&self) -> Vec<f64> {
        // let bbox = self.get_bbox();
        let deltax = self.bbox[2] - self.bbox[0];
        let deltay = self.bbox[3] - self.bbox[1];
        // let corners: Vec<f64> = Vec::new();
        let mut totals: Vec<f64> = vec![0., 0., 0., 0.];
        let mut ns: Vec<usize> = vec![0, 0, 0, 0];
        for i in (0..self.pts.len()).step_by(10) {
            let dx = (self.pts[i][0] - self.bbox[0]) / deltax;
            let dy = (self.pts[i][1] - self.bbox[1]) / deltay;
            if dx < 0.5 && dy < 0.5 {
                totals[0] += self.pts[i][2];
                ns[0] += 1;
            };
            if dx > 0.5 && dy < 0.5 {
                totals[1] += self.pts[i][2];
                ns[1] += 1;
            };
            if dx > 0.5 && dy > 0.5 {
                totals[2] += self.pts[i][2];
                ns[2] += 1;
            };
            if dx < 0.5 && dy > 0.5 {
                totals[3] += self.pts[i][2];
                ns[3] += 1;
            };
        }
        vec![
            totals[0] / ns[0] as f64,
            totals[1] / ns[1] as f64,
            totals[2] / ns[2] as f64,
            totals[3] / ns[3] as f64,
        ]
    }

    fn finalise(&self, impdigits: usize) -> io::Result<()> {
        if self.pts.is_empty() {
            io::stdout().write_all(&format!("{}\n", &self.cell).as_bytes())?;
            return Ok(());
        }
        info!("simplify {} points", self.pts.len());
        let bbox = &self.bbox;
        // let zavg = self.get_average_elevation();
        let cornerz = self.get_corner_elevations();
        let mut dt = startin::Triangulation::new();
        //-- insert 4 dummy corner TODO: better z-values would be nice
        let bufferbbox = 1000.0;
        let _ = dt.insert_one_pt(bbox[0] - bufferbbox, bbox[1] - bufferbbox, cornerz[0]);
        let _ = dt.insert_one_pt(bbox[2] + bufferbbox, bbox[1] - bufferbbox, cornerz[1]);
        let _ = dt.insert_one_pt(bbox[2] + bufferbbox, bbox[3] + bufferbbox, cornerz[2]);
        let _ = dt.insert_one_pt(bbox[0] - bufferbbox, bbox[3] + bufferbbox, cornerz[3]);
        //-- shuffle the input points
        let mut ids: Vec<usize> = (0..self.pts.len()).collect();
        ids.shuffle(&mut thread_rng());
        for i in ids {
            let p = &self.pts[i];
            if self.close_to_bbox(&p, 0.02) == false {
                let z2 = dt.interpolate_tin_linear(p[0], p[1]).unwrap();
                let e = (p[2] - z2).abs();
                if e > self.max_epsilon {
                    let _ = dt.insert_one_pt(p[0], p[1], p[2]);
                }
            }
        }
        //-- stream out the vertices
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
        info!(
            "{:.1?}% kept",
            (dt.number_of_vertices() - 4) as f64 / self.pts.len() as f64 * 100.
        );
        io::stdout().write_all(&format!("{}\n", &self.cell).as_bytes())?;
        Ok(())
    }
}

fn main() -> io::Result<()> {
    let cli = Cli::parse();
    env_logger::Builder::new()
        .filter_level(cli.verbose.log_level_filter())
        .init();

    // let pool = rayon::ThreadPoolBuilder::new()
    //     .num_threads(cli.cores)
    //     .build()
    //     .unwrap();
    // info!("Number of cores/threads used: {}", cli.cores);

    let mut surfaces: Vec<Surface> = vec![Surface::new(cli.vepsilon)];
    let mut impdigits: usize = usize::MAX;

    let stdin = std::io::stdin();
    for line in stdin.lock().lines() {
        let l = line.unwrap();
        if l.is_empty() {
            continue;
        }
        let s = surfaces.last_mut().unwrap();
        let ch = l.chars().next().unwrap();
        match ch {
            '#' => (),
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
                //-- find how many digits are stored in the first vertex
                //-- and use this for all
                if impdigits == usize::MAX {
                    let ls: Vec<&str> = l.split_whitespace().collect();
                    impdigits = ls[1].len() - 1 - ls[1].rfind('.').unwrap() as usize;
                    info!("impdigits={:?}", impdigits);
                }
                // surfaces.last_mut().unwrap().add_pt(parse_3_f64(&l));
                s.add_pt(parse_3_f64(&l));
            }
            'w' => {
                //-- sprinkled point
                io::stdout().write_all(&format!("{}\n", &l).as_bytes())?;
            }
            'x' => {
                //-- finalise a cell
                s.cell = l;
                info!("1. {:?}", surfaces.len());
                //-- wait till we have 10 in the vec, then use rayon
                if surfaces.len() >= 16 {
                    info!("trigger");
                    surfaces.par_iter().for_each(move |s| {
                        let _ = s.finalise(impdigits);
                    });
                    surfaces.clear();
                }
                surfaces.push(Surface::new(cli.vepsilon));
                info!("2. {:?}", surfaces.len());
            }
            _ => {
                error!("Wrongly formatted stream. Abort.");
                std::process::exit(1);
            }
        }
    }
    //-- write the ones not finalised
    warn!("size surfaces: {}", surfaces.len());
    surfaces.pop();
    surfaces.par_iter().for_each(move |s| {
        let _ = s.finalise(impdigits);
    });
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
