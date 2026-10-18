//-- ttsimpl

#[allow(dead_code)]
#[allow(unused_variables)]
#[macro_use]
extern crate log; //info/debug/error

extern crate startin;

use std::io::BufRead;
use std::io::{self, Write};

use rand::seq::SliceRandom;

use startin::interpolation;

use clap::Parser;
#[derive(Parser)]
#[command(name = "ttsimpl")]
#[command(about = "tintin -- simplify the terrain [ttsimpl]")]
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

    fn finalise(&mut self, out: &mut SpaOut) -> io::Result<()> {
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
        ids.shuffle(&mut rand::rng());
        for i in ids {
            let p = &self.pts[i];
            if self.close_to_bbox(&p, 0.02) == false {
                let locs = vec![[p[0], p[1]]];
                let mut re =
                    interpolation::interpolate(&interpolation::TIN {}, &mut dt, &locs);
                let z2 = re.remove(0).unwrap();
                let e = (p[2] - z2).abs();
                if e > self.max_epsilon {
                    let _ = dt.insert_one_pt(p[0], p[1], p[2]);
                }
            }
        }
        //-- stream out the vertices
        let allv = &dt.all_vertices();
        for i in 5..dt.number_of_vertices() {
            out.rec(&sst::stream::SpaRecord::V(
                allv[i][0], allv[i][1], allv[i][2],
            ))?;
        }
        info!(
            "{:.1?}% kept",
            (dt.number_of_vertices() - 4) as f64 / self.pts.len() as f64 * 100.
        );
        Ok(())
    }
}

//-- the `.spa` output of this tool: it mirrors the format of its input
//-- (a binary stream is simplified to a binary stream, an ASCII one to
//-- ASCII) and is buffered
struct SpaOut {
    w: std::io::BufWriter<std::io::StdoutLock<'static>>,
    binary: bool,
}

impl SpaOut {
    fn rec(&mut self, r: &sst::stream::SpaRecord) -> io::Result<()> {
        if self.binary == true {
            sst::stream::spa_write(&mut self.w, r)
        } else {
            self.w.write_all(r.to_line().as_bytes())
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        self.w.flush()
    }
}

fn main() -> io::Result<()> {
    let cli = Cli::parse();
    env_logger::Builder::new()
        .filter_level(cli.verbose.log_level_filter())
        .init();

    let mut s = Surface::new(cli.vepsilon);

    //-- a binary input is answered with a binary output; the records are
    //-- the same either way
    let stdin = std::io::stdin();
    let mut r = std::io::BufReader::new(stdin.lock());
    let binary = sst::stream::peek_is_binary(&mut r, sst::stream::SPA_MAGIC)?;
    if binary == true {
        sst::stream::skip_magic(&mut r, sst::stream::SPA_MAGIC)?;
    }
    let stdout = io::stdout();
    let mut out = SpaOut {
        w: std::io::BufWriter::with_capacity(1 << 20, stdout.lock()),
        binary,
    };
    if binary == true {
        sst::stream::spa_write_magic(&mut out.w)?;
    }

    if binary == true {
        while let Some(rec) = sst::stream::spa_read(&mut r)? {
            handle_record(&mut s, &mut out, &rec)?;
        }
    } else {
        for line in r.lines() {
            let l = line.unwrap();
            if l.is_empty() {
                continue;
            }
            let ch = l.chars().next().unwrap();
            if ch == '#' {
                continue;
            }
            match parse_sp_record(&l) {
                Some(rec) => handle_record(&mut s, &mut out, &rec)?,
                None => {
                    error!("Wrongly formatted stream. Abort.");
                    std::process::exit(1);
                }
            }
        }
    }
    out.flush()?;

    info!("✅");
    Ok(())
}

/// One `.spa` record, from an ASCII line (the binary path decodes the same
/// enum directly).
fn handle_record(s: &mut Surface, out: &mut SpaOut, rec: &sst::stream::SpaRecord) -> io::Result<()> {
    use sst::stream::SpaRecord;
    match rec {
        SpaRecord::N(_) | SpaRecord::C(_) | SpaRecord::S(_) | SpaRecord::B(_) => {
            out.rec(rec)?;
        }
        SpaRecord::W(..) => {
            //-- sprinkled points are not simplified, just passed on
            out.rec(rec)?;
        }
        SpaRecord::V(x, y, z) => {
            s.add_pt(vec![*x, *y, *z]);
        }
        SpaRecord::X(_, _) => {
            //-- finalise a cell: simplify what was buffered, then pass it on
            if s.is_empty() == false {
                s.finalise(out)?;
                s.clear();
            }
            out.rec(rec)?;
        }
    }
    Ok(())
}

/// Parse one ASCII `.spa` line into a record; `None` when malformed.
fn parse_sp_record(l: &str) -> Option<sst::stream::SpaRecord> {
    use sst::stream::SpaRecord;
    let mut it = l.split_whitespace();
    let tag = it.next()?;
    Some(match tag {
        "n" => SpaRecord::N(it.next()?.parse::<u64>().ok()?),
        "c" => SpaRecord::C(it.next()?.parse::<u32>().ok()?),
        "s" => SpaRecord::S(it.next()?.parse::<u64>().ok()?),
        "b" => SpaRecord::B([
            it.next()?.parse::<f64>().ok()?,
            it.next()?.parse::<f64>().ok()?,
            it.next()?.parse::<f64>().ok()?,
            it.next()?.parse::<f64>().ok()?,
        ]),
        "v" => SpaRecord::V(
            it.next()?.parse::<f64>().ok()?,
            it.next()?.parse::<f64>().ok()?,
            it.next()?.parse::<f64>().ok()?,
        ),
        "w" => SpaRecord::W(
            it.next()?.parse::<f64>().ok()?,
            it.next()?.parse::<f64>().ok()?,
            it.next()?.parse::<f64>().ok()?,
        ),
        "x" => SpaRecord::X(
            it.next()?.parse::<u32>().ok()?,
            it.next()?.parse::<u32>().ok()?,
        ),
        _ => return None,
    })
}
