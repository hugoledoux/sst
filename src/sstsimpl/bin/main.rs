//-- sstsimpl

#[allow(dead_code)]
#[allow(unused_variables)]
#[macro_use]
extern crate log; //info/debug/error

extern crate startin;

use std::io::BufRead;
use std::io::{self, Write};

use rand::seq::SliceRandom;
use rand::thread_rng;

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

struct Surface {
    pts: Vec<Vec<f64>>,
    max_epsilon: f64,
    bbox: Vec<f64>,
}

impl Surface {
    fn new(epsilon: f64) -> Surface {
        let l = Vec::new();
        let b = vec![f64::MAX, f64::MAX, f64::MIN, f64::MIN];
        Surface {
            pts: l,
            max_epsilon: epsilon,
            bbox: b,
        }
    }

    fn clear(&mut self) {
        self.pts.clear();
        self.bbox = vec![f64::MAX, f64::MAX, f64::MIN, f64::MIN];
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

    // fn get_4_extremes(&self) -> Vec<usize> {
    //     let mut left = 0;
    //     let mut bottom = 0;
    //     let mut right = 0;
    //     let mut top = 0;
    //     for i in 1..self.pts.len() {
    //         if self.pts[i][0] < self.pts[left][0] {
    //             left = i;
    //         }
    //         if self.pts[i][0] > self.pts[right][0] {
    //             right = i;
    //         }
    //         if self.pts[i][1] < self.pts[bottom][1] {
    //             bottom = i;
    //         }
    //         if self.pts[i][1] > self.pts[top][1] {
    //             top = i;
    //         }
    //     }
    //     vec![left, bottom, right, top]
    // }

    // fn get_average_elevation(&self) -> f64 {
    //     let mut total = 0_f64;
    //     let mut n: usize = 0;
    //     for i in (0..self.pts.len()).step_by(10) {
    //         total += self.pts[i][2] as f64;
    //         n += 1;
    //     }
    //     total / n as f64
    // }

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

    fn finalise(&mut self) -> io::Result<()> {
        // info!("->: {}", self.pts.len());
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
        let impdigits = 3; // TODO : impdigit to store in the stream?
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
            'w' => {
                //-- sprinkled point
                io::stdout().write_all(&format!("{}\n", &l).as_bytes())?;
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
