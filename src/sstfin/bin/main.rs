use std::fmt;
use std::fs::File;
use std::io::Seek;
use std::io::{BufRead, BufReader};

#[derive(Clone)]
pub struct Point {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({}, {}, {})", self.x, self.y, self.z)
    }
}

fn main() {
    // stream1();

    let f = File::open("/Users/hugo/code/dt-comparisons/data/5.txt").expect("Unable to open file");

    //-- pass #1
    let mut bbox = pass_1(&f);
    // println!("bbox {:?}", bbox);
    let deltax = bbox[2] - bbox[0];
    let deltay = bbox[3] - bbox[1];
    if deltax >= deltay {
        bbox[3] = bbox[1] + deltax;
    } else {
        bbox[2] = bbox[0] + deltay;
    }
    println!("bbox {:?}", bbox);

    //-- pass #2
    let cellsize: usize = 10;
    let mut g: Vec<Vec<usize>> = pass_2(&f, &bbox, cellsize);

    //-- pass #3
    pass_3(&f, &bbox, cellsize, &mut g);
}

fn pass_3(mut f: &File, bbox: &Vec<f64>, cellsize: usize, g: &mut Vec<Vec<usize>>) {
    let width: usize = ((bbox[2] - bbox[0]) / cellsize as f64).ceil() as usize;
    let height: usize = ((bbox[3] - bbox[1]) / cellsize as f64).ceil() as usize;
    let mut gpts: Vec<Vec<Vec<Point>>> = vec![vec![Vec::new(); height]; width];
    let _re = f.seek(std::io::SeekFrom::Start(0)); //-- reset to begining of the file
    let f = BufReader::new(f);
    for l in f.lines() {
        let l = l.expect("Unable to read line");
        let v: Vec<f64> = l.split(' ').map(|s| s.parse().unwrap()).collect();
        let gxy: (usize, usize) = get_gx_gy(v[0], v[1], bbox[0], bbox[1], cellsize);
        println!("{}--{}", gxy.0, gxy.1);
        g[gxy.0][gxy.1] -= 1;
        gpts[gxy.0][gxy.1].push(Point {
            x: v[0],
            y: v[1],
            z: v[2],
        });
        if g[gxy.0][gxy.1] == 0 {
            println!("FINALISATION OF CELL {}--{}", gxy.0, gxy.1);
        }
    }
}

// fn finalise_cell(lspts: &mut <Vec<Point>>)
// {

// }

fn pass_2(mut f: &File, bbox: &Vec<f64>, cellsize: usize) -> Vec<Vec<usize>> {
    let width: usize = ((bbox[2] - bbox[0]) / cellsize as f64).ceil() as usize;
    let height: usize = ((bbox[3] - bbox[1]) / cellsize as f64).ceil() as usize;
    let mut g: Vec<Vec<usize>> = vec![vec![0; height]; width];
    let _re = f.seek(std::io::SeekFrom::Start(0)); //-- reset to begining of the file
    let f = BufReader::new(f);
    for l in f.lines() {
        let l = l.expect("Unable to read line");
        let v: Vec<f64> = l.split(' ').map(|s| s.parse().unwrap()).collect();
        let gxy: (usize, usize) = get_gx_gy(v[0], v[1], bbox[0], bbox[1], cellsize);
        // println!("{}--{}", gxy.0, gxy.1);
        g[gxy.0][gxy.1] += 1;
    }
    g
}

fn get_gx_gy(x: f64, y: f64, minx: f64, miny: f64, cellsize: usize) -> (usize, usize) {
    (
        ((x - minx) / cellsize as f64) as usize,
        ((y - miny) / cellsize as f64) as usize,
    )
}

fn pass_1(f: &File) -> Vec<f64> {
    let f = BufReader::new(f);
    let mut bbox: Vec<f64> = Vec::new();
    bbox.push(std::f64::MAX);
    bbox.push(std::f64::MAX);
    bbox.push(std::f64::MIN);
    bbox.push(std::f64::MIN);
    for l in f.lines() {
        let l = l.expect("Unable to read line");
        let v: Vec<f64> = l.split(' ').map(|s| s.parse().unwrap()).collect();
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
    bbox
}

// fn stream1() {
//     let stdin = std::io::stdin();
//     for line in stdin.lock().lines() {
//         println!("{}", line.unwrap());
//     }
// }
