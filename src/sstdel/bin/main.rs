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
    let mut _totalpts: usize = 0;
    // let mut cellsize: usize = 0;
    // let mut bbox: [f64; 2] = [std::f64::MIN, std::f64::MIN];
    // let mut gwidth: usize = 0;
    // let mut gheight: usize = 0;

    info!("Init DT");
    let mut dt = startin::Triangulation::new();

    //PUTBACK let stdin = std::io::stdin();
    //PUTBACK for line in stdin.lock().lines() {
    let fi = File::open("/Users/hugo/projects/sst/data/rect1.stream").expect("Unable to open file");
    // File::open("/Users/hugo/projects/sst/data/square400.stream").expect("Unable to open file");
    let f = BufReader::new(fi);
    let mut count: usize = 0;
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
                _totalpts = l
                    .split_whitespace()
                    .last()
                    .unwrap()
                    .parse::<usize>()
                    .unwrap()
            }
            'r' => {
                //-- cellsize
                let r = l
                    .split_whitespace()
                    .last()
                    .unwrap()
                    .parse::<usize>()
                    .unwrap();
                dt.set_cellsize(r);
            }
            'd' => {
                //-- dimension grid
                let re = parse_2_usize(&l);
                dt.set_grid_dimensions(re.0, re.1);
            }
            'b' => {
                //-- bbox
                let re = parse_2_f64(&l);
                dt.set_bbox(re.0, re.1);
                info!("Writing GeoJSON file to disk: /Users/hugo/temp/sstout/z.grid.geojson");
                let _re =
                    dt.write_geojson_grid("/Users/hugo/temp/sstout/z.grid.geojson".to_string());
            }
            'v' => {
                //-- vertex
                println!("{}", l);
                count += 1;
                let v = parse_3_f64(&l);
                let _re = dt.insert_one_pt_with_grid(v.0, v.1, v.2);
            }
            'x' => {
                println!("{}", l);
                //-- finalise a cell
                let re = parse_2_usize(&l);
                let _re = dt.finalise_cell(re.0, re.1);
                let fout = format!("/Users/hugo/temp/sstout/c-{}-{}.geojson", re.0, re.1);
                let _re = dt.write_geojson_triangles(fout.to_string());
                if ((re.0 == 3) && (re.1 == 1)) {
                    println!("yo");
                }
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

    // info!("Writing GeoJSON file to disk: /Users/hugo/temp/z.grid.geojson");
    // let _re = dt.write_geojson_grid("/Users/hugo/temp/z.grid.geojson".to_string());

    info!("Writing GeoJSON file to disk: /Users/hugo/temp/sstout/z.geojson");
    let _re = dt.write_geojson_triangles("/Users/hugo/temp/sstout/z.geojson".to_string());

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
