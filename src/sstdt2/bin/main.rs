//-- sstdt2

#[allow(dead_code)]
#[allow(unused_variables)]
mod startin;

#[macro_use]
extern crate log; //info/debug/error

use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::{self, Write};

use clap::App;

fn main() -> io::Result<()> {
    let _matches = App::new("sstdt2")
        .version("0.1")
        .about("streaming startin -- Delaunay triangulation")
        .arg("--stars...       'output stars instead of .sma'")
        .arg("--both...        'output both vertices and stars'")
        .get_matches();

    env_logger::init();

    let mut _totalpts: usize = 0;

    info!("Init DT");
    let mut dt = startin::Triangulation::new();

    // if matches.occurrences_of("stars") > 0 {
    //     dt.set_outputmode(Outputmode::Stars);
    // }

    // if matches.occurrences_of("both") > 0 {
    //     dt.set_outputmode(Outputmode::Both);
    // }

    //----- reading from stdin -----//
    // let stdin = std::io::stdin();
    // for line in stdin.lock().lines() {
    //     let l = line.unwrap();
    //----- reading from stdin -----//

    //----- reading from file -----//
    // let fi = File::open("/Users/hugo/projects/sst/data/s400_50.spa").expect("Unable to open file");
    let fi = File::open("/Users/hugo/projects/sst/data/rect1_50.spa").expect("Unable to open file");
    // let fi = File::open("/Users/hugo/projects/sst/data/test1.spa").expect("Unable to open file");
    let f = BufReader::new(fi);
    for l in f.lines() {
        let l = l.expect("Unable to read line");
        //----- reading from file -----//

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
            's' => {
                //-- cellsize
                let s = l
                    .split_whitespace()
                    .last()
                    .unwrap()
                    .parse::<usize>()
                    .unwrap();
                dt.set_cellsize(s);
            }
            'c' => {
                //-- dimension grid (always square cXc)
                let c = l
                    .split_whitespace()
                    .last()
                    .unwrap()
                    .parse::<usize>()
                    .unwrap();
                dt.set_grid_dimensions(c);
            }
            'b' => {
                //-- bbox
                let re = parse_4_f64(&l);
                dt.set_bbox(re.0, re.1, re.2, re.3);
                io::stdout().write_all(
                    &format!("b {:.3} {:.3} {:.3} {:.3}\n", re.0, re.1, re.2, re.3).as_bytes(),
                )?;
            }
            'v' => {
                //-- vertex
                // println!("=>{}", l);
                // count += 1;
                let v = parse_3_f64(&l);
                let _re = dt.insert_one_pt_with_grid(v.0, v.1, v.2);
                // if v.0 > 195.0 && v.0 < 196.0 && v.1 > 2.0 && v.1 < 3.0 {
                //     //     // println!("that one");
                //     let fout = format!("/Users/hugo/temp/sstout/tmp.geojson");
                //     let _re = dt.write_geojson_triangles(fout.to_string());
                //     //     dt.print_ds();
                // }
            }
            'x' => {
                //-- finalise a cell
                // println!("=>{}", l);
                let re = parse_2_usize(&l);
                if (re.0 == 1) && (re.1 == 0) {
                    let fout = format!("/Users/hugo/temp/sstout/c-{}-{}-b.geojson", re.0, re.1);
                    let _re = dt.write_geojson_triangles(fout.to_string());
                }

                let _re = dt.finalise_qtcell(re.0, re.1);
                if (re.0 == 1) && (re.1 == 0) {
                    let fout = format!("/Users/hugo/temp/sstout/c-{}-{}.geojson", re.0, re.1);
                    let _re = dt.write_geojson_triangles(fout.to_string());
                }
            }
            _ => {
                error!("Wrongly formatted stream. Abort.");
                std::process::exit(1);
            }
        }
    }

    // println!("==> {:?}", dt.incident_triangles_to_vertex(3));

    // info!("Writing GeoJSON file to disk: /Users/hugo/temp/sstout/z.grid.geojson");
    // let _re = dt.write_geojson_grid("/Users/hugo/temp/sstout/z.grid.geojson".to_string());

    // let _re = dt.write_geojson_triangles("/Users/hugo/temp/sstout/tr.geojson".to_string());

    // dt.print_ds();

    info!("Finished reading the stream");
    info!("dt.number_of_vertices() = {}", dt.number_of_vertices());

    // info!("max # points in DT during process: {}", dt.max);
    // let _x = dt.finalise_leftover_triangles();

    info!("✅");
    Ok(())
}

fn parse_2_usize(l: &String) -> (usize, usize) {
    let ls: Vec<&str> = l.split_whitespace().collect();
    let a: usize = ls[1].parse::<usize>().unwrap();
    let b: usize = ls[2].parse::<usize>().unwrap();
    (a, b)
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
