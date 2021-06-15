//-- sstval

#[allow(dead_code)]
#[allow(unused_variables)]
#[macro_use]
extern crate log; //info/debug/error

use hashbrown::HashMap;
use std::io::BufRead;
use std::io::{self};

use clap::App;

fn main() -> io::Result<()> {
    let _matches = App::new("sstval")
        .version("0.2")
        .about("validate the sobj file")
        .get_matches();

    env_logger::init();

    let mut total_vs: usize = 0;
    let mut total_ts: usize = 0;

    let mut d: HashMap<usize, (f64, f64, f64)> = HashMap::new();

    let stdin = std::io::stdin();
    for line in stdin.lock().lines() {
        let l = line.unwrap();
        if l.is_empty() {
            continue;
        }
        let ch = l.chars().next().unwrap();
        match ch {
            '#' => continue,
            'b' => {
                let bbox = parse_4_f64(&l);
                info!("bbox: {:?}", bbox);
            }
            'v' => {
                //-- a vertex
                let v = parse_3_f64(&l);
                let _re = d.insert(total_vs, v);
                total_vs += 1;
            }
            'f' => {
                //-- a triangle
                // let re = parse_3_usize(&l);
                total_ts += 1;
            }
            'x' => {
                //-- finalise a vertex
                let ls: Vec<&str> = l.split_whitespace().collect();
                let a: usize = ls[1].parse::<usize>().unwrap();
                d.remove(&a);
            }
            _ => {
                error!("Wrongly formatted stream. Abort.");
                std::process::exit(1);
            }
        }
    }
    info!("number of vertices: {}", total_vs);
    info!("number of triangles: {}", total_ts);
    info!("✅");
    Ok(())
}

fn parse_3_f64(l: &String) -> (f64, f64, f64) {
    let ls: Vec<&str> = l.split_whitespace().collect();
    let a: f64 = ls[1].parse::<f64>().unwrap();
    let b: f64 = ls[2].parse::<f64>().unwrap();
    let c: f64 = ls[3].parse::<f64>().unwrap();
    (a, b, c)
}

fn parse_4_f64(l: &String) -> (f64, f64, f64, f64) {
    let ls: Vec<&str> = l.split_whitespace().collect();
    let a: f64 = ls[1].parse::<f64>().unwrap();
    let b: f64 = ls[2].parse::<f64>().unwrap();
    let c: f64 = ls[3].parse::<f64>().unwrap();
    let d: f64 = ls[4].parse::<f64>().unwrap();
    (a, b, c, d)
}
