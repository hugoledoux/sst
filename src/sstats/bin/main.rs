// sstats

use std::fmt;
use std::fs::File;
use std::path::Path;

use std::io::Seek;
use std::io::{BufRead, BufReader};

use std::collections::HashMap;

extern crate las;

use las::Read;

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

    let arg = std::env::args().skip(1).next();
    let s = match arg {
        Some(p) => p,
        None => {
            eprintln!("Must provide a path to a lidar file.");
            std::process::exit(1);
        }
    };

    let path = Path::new(&s);
    if path.extension().unwrap() == "las"
        || path.extension().unwrap() == "laz"
        || path.extension().unwrap() == "LAS"
        || path.extension().unwrap() == "LAZ"
    {
        input_las(&path);
    } else if path.extension().unwrap() == "txt" || path.extension().unwrap() == "xyz" {
        input_csv(&path);
    }
}

fn input_las(path: &Path) {
    let reader = las::Reader::from_path(path.clone()).expect("Unable to open reader");
    let totalpts = reader.header().number_of_points();
    let b = reader.header().bounds();
    let bbox = vec![b.min.x, b.min.y, b.max.x, b.max.y];
    println!("count:   {:>7}", totalpts);
    println!("x-extent: {:>5.3}", bbox[2] - bbox[0]);
    println!("y-extent: {:>5.3}", bbox[3] - bbox[1]);
    //-- make it a square to have a quadtree
    let size = if (bbox[3] - bbox[1]) > (bbox[2] - bbox[0]) {
        bbox[3] - bbox[1]
    } else {
        bbox[2] - bbox[0]
    };
    println!("| cells | resolution | avg/cell | max/cell |");
    let re = stats_per_cell_las(&path, &bbox, size);
    for each in re {
        println!(
            " {:3}x{:<3}    {:3}m        {:6}     {:6}",
            each.0, each.0, each.1, each.2 as u32, each.3
        );
    }
}

fn input_csv(path: &Path) {
    let f = File::open(path).expect("Unable to open file");
    let (bbox, totalpts) = get_count_extent(&f);
    println!("count:   {:>7}", totalpts);
    println!("x-extent: {:>5.3}", bbox[2] - bbox[0]);
    println!("y-extent: {:>5.3}", bbox[3] - bbox[1]);
    //-- make it a square to have a quadtree
    let size = if (bbox[3] - bbox[1]) > (bbox[2] - bbox[0]) {
        bbox[3] - bbox[1]
    } else {
        bbox[2] - bbox[0]
    };
    println!("| cells | resolution | avg/cell | max/cell |");
    let re = stats_per_cell_csv(&f, &bbox, size);
    for each in re {
        println!(
            " {:3}x{:<3}    {:3}m        {:6}     {:6}",
            each.0, each.0, each.1, each.2 as u32, each.3
        );
    }
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

fn stats_per_cell_csv(mut f: &File, bbox: &Vec<f64>, size: f64) -> Vec<(usize, usize, f64, usize)> {
    // (cells | resolution | avg/cell | max/cell)
    let depths: Vec<usize> = vec![3, 4, 5, 6, 7, 8, 9];
    let mut grids = HashMap::new();
    let mut results = Vec::new();
    for depth in &depths {
        let cellno: usize = (2_u32).pow(*depth as u32) as usize;
        grids.insert(depth, vec![vec![0; cellno]; cellno]);
    }
    for depth in &depths {
        let _re = f.seek(std::io::SeekFrom::Start(0)); //-- reset to begining of the file
        let f = BufReader::new(f);

        let cellno: usize = (2_u32).pow(*depth as u32) as usize;
        let resolution = (size / cellno as f64).ceil() as usize;
        for l in f.lines() {
            let l = l.expect("Unable to read line");
            let v: Vec<f64> = l.split(' ').map(|s| s.parse().unwrap()).collect();
            let mut gxy: (usize, usize) = get_gx_gy(v[0], v[1], bbox[0], bbox[1], resolution);
            if gxy.0 == cellno {
                gxy.0 = cellno - 1;
            }
            if gxy.1 == cellno {
                gxy.1 = cellno - 1;
            }
            grids.get_mut(&depth).unwrap()[gxy.0][gxy.1] += 1;
        }
        let mut max: usize = 0;
        let mut count: usize = 0;
        let mut nonempty: usize = 0;
        for i in 0..cellno {
            for j in 0..cellno {
                if grids[&depth][i][j] > 0 {
                    nonempty += 1;
                    count += grids[&depth][i][j];
                }
                if grids[&depth][i][j] > max {
                    max = grids[&depth][i][j];
                }
            }
        }
        let avg: f64 = count as f64 / nonempty as f64;
        results.push((cellno, resolution, avg, max));
    }
    results
}

fn stats_per_cell_las(path: &Path, bbox: &Vec<f64>, size: f64) -> Vec<(usize, usize, f64, usize)> {
    // (cells | resolution | avg/cell | max/cell)
    let depths: Vec<usize> = vec![3, 4, 5, 6, 7, 8, 9];
    let mut grids = HashMap::new();
    let mut results = Vec::new();
    for depth in &depths {
        let cellno: usize = (2_u32).pow(*depth as u32) as usize;
        grids.insert(depth, vec![vec![0; cellno]; cellno]);
    }
    for depth in &depths {
        let mut reader = las::Reader::from_path(path).expect("Unable to open reader");
        let cellno: usize = (2_u32).pow(*depth as u32) as usize;
        let resolution = (size / cellno as f64).ceil() as usize;
        for each in reader.points() {
            let p = each.unwrap();
            let mut gxy: (usize, usize) = get_gx_gy(p.x, p.y, bbox[0], bbox[1], resolution);
            if gxy.0 == cellno {
                gxy.0 = cellno - 1;
            }
            if gxy.1 == cellno {
                gxy.1 = cellno - 1;
            }
            grids.get_mut(&depth).unwrap()[gxy.0][gxy.1] += 1;
        }
        let mut max: usize = 0;
        let mut count: usize = 0;
        let mut nonempty: usize = 0;
        for i in 0..cellno {
            for j in 0..cellno {
                if grids[&depth][i][j] > 0 {
                    nonempty += 1;
                    count += grids[&depth][i][j];
                }
                if grids[&depth][i][j] > max {
                    max = grids[&depth][i][j];
                }
            }
        }
        let avg: f64 = count as f64 / nonempty as f64;
        results.push((cellno, resolution, avg, max));
    }
    results
}

fn get_gx_gy(x: f64, y: f64, minx: f64, miny: f64, resolution: usize) -> (usize, usize) {
    (
        ((x - minx) / resolution as f64) as usize,
        ((y - miny) / resolution as f64) as usize,
    )
}
