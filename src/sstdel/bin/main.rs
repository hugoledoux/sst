//-- sstdel

#[allow(dead_code)]
#[allow(unused_variables)]
mod startin;

#[macro_use]
extern crate log; //info/debug/error

use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader};

fn main() -> io::Result<()> {
    env_logger::init();
    let mut totalpts: usize = 0;
    let mut cellsize: usize = 0;
    let mut bbox: [f64; 2] = [std::f64::MIN, std::f64::MIN];
    let mut gwidth: usize = 0;
    let mut gheight: usize = 0;

    let mut dt = startin::Triangulation::new(0, 0, 0, std::f64::MAX, std::f64::MAX);
    // let mut dt: startin::Triangulation;

    info!("Init DT");

    //PUTBACK let stdin = std::io::stdin();
    //PUTBACK for line in stdin.lock().lines() {
    let fi =
        File::open("/Users/hugo/projects/sst/data/square400.stream").expect("Unable to open file");
    let f = BufReader::new(fi);
    for l in f.lines() {
        let l = l.expect("Unable to read line");
        //PUTBACK let l = line.unwrap();
        // println!("=> {}", l);
        if l.is_empty() {
            continue;
        }
        let ch = l.chars().next().unwrap();
        match ch {
            '#' => continue,
            'n' => {
                //-- number of points
                totalpts = l
                    .split_whitespace()
                    .last()
                    .unwrap()
                    .parse::<usize>()
                    .unwrap()
            }
            'c' => {
                //-- cellsize
                cellsize = l
                    .split_whitespace()
                    .last()
                    .unwrap()
                    .parse::<usize>()
                    .unwrap()
            }
            'd' => {
                //-- dimension grid
                let re = parse_2_usize(&l);
                gwidth = re.0;
                gheight = re.1;
            }
            'b' => {
                //-- bbox
                let re = parse_2_f64(&l);
                bbox[0] = re.0;
                bbox[1] = re.1;
                //-- init the DT
                dt = startin::Triangulation::new(gwidth, gheight, cellsize, re.0, re.1);
            }
            'v' => {
                //-- vertex
                let v = parse_3_f64(&l);
                let _re = dt.insert_one_pt(v.0, v.1, v.2);
            }
            'x' => {
                //-- finalise a cell
                let re = parse_2_usize(&l);
                info!("Cell {}--{} finalised", re.0, re.1,);
                let _re = dt.finalise_cell(re.0, re.1);
            }
            _ => {
                error!("Wrongly formatted stream. Abort.");
                std::process::exit(1);
            }
        }
    }
    info!("Finished reading the stream");
    info!("dt.number_of_vertices() = {}", dt.number_of_vertices());

    // println!("{}", dt.printme(false));
    // std::process::exit(1);

    let _x = dt.finalise_leftover_triangles();
    Ok(())
}

fn parse_2_usize(l: &String) -> (usize, usize) {
    let ls: Vec<&str> = l.split_whitespace().collect();
    let a: usize = ls[1].parse::<usize>().unwrap();
    let b: usize = ls[2].parse::<usize>().unwrap();
    (a, b)
}

fn parse_2_f64(l: &String) -> (f64, f64) {
    let ls: Vec<&str> = l.split_whitespace().collect();
    let a: f64 = ls[1].parse::<f64>().unwrap();
    let b: f64 = ls[2].parse::<f64>().unwrap();
    (a, b)
}

fn parse_3_f64(l: &String) -> (f64, f64, f64) {
    let ls: Vec<&str> = l.split_whitespace().collect();
    let a: f64 = ls[1].parse::<f64>().unwrap();
    let b: f64 = ls[2].parse::<f64>().unwrap();
    let c: f64 = ls[3].parse::<f64>().unwrap();
    (a, b, c)
}
