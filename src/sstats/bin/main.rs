use std::env;
use std::fmt;
use std::fs::File;

use std::io::Seek;
use std::io::{BufRead, BufReader};

#[macro_use]
extern crate log; //info/debug/error

#[derive(Clone)]
pub struct Point {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point {
    pub fn new() -> Point {
        Point {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }
}

impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({}, {}, {})", self.x, self.y, self.z)
    }
}

fn main() {
    env_logger::init();

    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        // error!("Input parameter not given: cellsize. Abort.");
        // eprintln!("Input parameter not given: cellsize");
        std::process::exit(1);
    }
    let f = File::open(&args[1]).expect("Unable to open file");

    let (bbox, totalpts) = get_count_extent(&f);
    // info!("bbox={:?}", bbox);
    println!("count:   {:>8}", totalpts);
    println!("x-extent: {:>5.3}", bbox[2] - bbox[0]);
    println!("y-extent: {:>5.3}", bbox[3] - bbox[1]);

    //-- make it a square to have a quadtree
    let size = if (bbox[3] - bbox[1]) > (bbox[2] - bbox[0]) {
        (bbox[3] - bbox[1])
    } else {
        (bbox[2] - bbox[0])
    };

    println!("| cells | resolution | avg/cell | max/cell |");
    for i in 3..10 {
        let cellno: usize = (2_u32).pow(i) as usize;
        let resolution = (size / cellno as f64).ceil() as usize;
        let (avg, max) = stats_per_cell(&f, &bbox, cellno, resolution);
        println!(
            " {:3}x{:<3}    {:3}m        {:6}     {:6}",
            cellno, cellno, resolution, avg as u32, max
        );
        // println!("avg {}", avg);
        // println!("max {}", max);
    }

    // //-- pass #2
    // let mut g: Vec<Vec<usize>> = pass_2(&f, &bbox, cellsize, &mut sprinkled);
    // info!("Second pass ✅");
}

fn get_count_extent(f: &File) -> (Vec<f64>, usize) {
    let f = BufReader::new(f);
    let mut bbox: Vec<f64> = Vec::new();
    bbox.push(std::f64::MAX);
    bbox.push(std::f64::MAX);
    bbox.push(std::f64::MIN);
    bbox.push(std::f64::MIN);
    let mut n: usize = 0;
    for l in f.lines() {
        let l = l.expect("Unable to read line");
        let v: Vec<f64> = l.split(' ').map(|s| s.parse().unwrap()).collect();
        n += 1;
        //-- minxy
        for i in 0..2 {
            if v[i] < bbox[i] {
                bbox[i] = v[i];
            }
        }
        //-- maxxy
        for i in 0..2 {
            if v[i] > bbox[i + 2] {
                bbox[i + 2] = v[i];
            }
        }
    }
    (bbox, n)
}

fn stats_per_cell(mut f: &File, bbox: &Vec<f64>, cellno: usize, resolution: usize) -> (f64, usize) {
    let mut g: Vec<Vec<usize>> = vec![vec![0; cellno]; cellno];
    let _re = f.seek(std::io::SeekFrom::Start(0)); //-- reset to begining of the file
    let f = BufReader::new(f);
    for l in f.lines() {
        let l = l.expect("Unable to read line");
        let v: Vec<f64> = l.split(' ').map(|s| s.parse().unwrap()).collect();
        let mut gxy: (usize, usize) = get_gx_gy(v[0], v[1], bbox[0], bbox[1], resolution);
        // println!("{}--{}", gxy.0, gxy.1);
        if gxy.0 == cellno {
            gxy.0 = cellno - 1;
        }
        if gxy.1 == cellno {
            gxy.1 = cellno - 1;
        }
        g[gxy.0][gxy.1] += 1;
    }
    let mut max: usize = 0;
    let mut count: usize = 0;
    let mut nonempty: usize = 0;
    for i in 0..cellno {
        for j in 0..cellno {
            if g[i][j] > 0 {
                nonempty += 1;
                count += g[i][j];
            }
            if g[i][j] > max {
                max = g[i][j];
            }
        }
    }
    let mut avg: f64 = count as f64 / nonempty as f64;
    (avg, max)
}

fn get_gx_gy(x: f64, y: f64, minx: f64, miny: f64, resolution: usize) -> (usize, usize) {
    (
        ((x - minx) / resolution as f64) as usize,
        ((y - miny) / resolution as f64) as usize,
    )
}
