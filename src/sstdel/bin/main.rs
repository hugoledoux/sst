//-- sstdel

#[allow(dead_code)]
#[allow(unused_variables)]
mod startin;

#[macro_use]
extern crate log; //info/debug/error
use std::collections::HashSet;

use std::io::BufRead;
use std::io::{self, Write};

fn main() -> io::Result<()> {
    env_logger::init();
    let mut totalpts: usize = 0;
    let mut cellsize: usize = 0;
    let mut bbox: [f64; 2] = [std::f64::MIN, std::f64::MIN];
    let mut gpts: Vec<Vec<HashSet<usize>>> = Vec::new();

    let mut dt = startin::Triangulation::new();
    info!("Init DT");

    let stdin = std::io::stdin();
    for line in stdin.lock().lines() {
        let l = line.unwrap();
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
                gpts.resize(re.0, Vec::new());
                for each in &mut gpts {
                    each.resize(re.1, HashSet::new());
                }
            }
            'b' => {
                //-- bbox
                let re = parse_2_f64(&l);
                bbox[0] = re.0;
                bbox[1] = re.1;
            }
            'v' => {
                //-- vertex
                let v = parse_3_f64(&l);
                let re = dt.insert_one_pt(v.0, v.1, v.2);
                match re {
                    Ok(x) => {
                        let g = get_gx_gy(v.0, v.1, bbox[0], bbox[1], cellsize);
                        gpts[g.0][g.1].insert(x);
                    }
                    Err(_x) => continue,
                };
            }
            'x' => {
                //-- finalise a cell
                let re = parse_2_usize(&l);
                if gpts[re.0][re.1].is_empty() == false {
                    let mut gbbox: [f64; 4] = [0.0, 0.0, 0.0, 0.0];
                    get_cell_bbox(re.0, re.1, &bbox, cellsize, &mut gbbox);
                    // println!("re={}-{}", re.0, re.1);
                    // println!("box={:?}", bbox);
                    // println!("cellsize={}", cellsize);
                    // println!("gbbox={:?}", gbbox);
                    finalise_cell(&mut dt, &mut gpts, (re.0, re.1), &gbbox);
                }
            }
            _ => {
                error!("Wrongly formatted stream. Abort.");
                std::process::exit(1);
            }
        }
    }
    info!("Finished reading the stream");
    let mut total: usize = 0;
    for w in &gpts {
        for h in w {
            total += h.len();
        }
    }
    info!("Writing the {} vertices left in the DT", total);
    info!("DT # points: {}", dt.number_of_vertices());
    //-- write the leftovers
    for w in &gpts {
        for h in w {
            for each in h.iter() {
                let p = dt.get_point(*each).unwrap();
                io::stdout().write_all(
                    &format!(
                        "v {} {} {} {} {:?}\n",
                        *each,
                        p[0],
                        p[1],
                        p[2],
                        dt.adjacent_vertices_to_vertex(*each).unwrap()
                    )
                    .as_bytes(),
                )?;
            }
        }
    }
    Ok(())
}

fn finalise_cell(
    dt: &mut startin::Triangulation,
    gpts: &mut Vec<Vec<HashSet<usize>>>,
    cellid: (usize, usize),
    gbbox: &[f64],
) -> io::Result<()> {
    let cell = &mut gpts[cellid.0][cellid.1];
    info!(
        "Cell {}--{} finalised ({} vertices)",
        cellid.0,
        cellid.1,
        cell.len()
    );
    //-- verify which vertices can be finalised/flushed
    let mut finpts: HashSet<usize> = HashSet::new();
    for theid in cell.iter() {
        // let p = dt.get_point(each);
        let re = dt.adjacent_vertices_to_vertex(*theid).unwrap();
        let mut fin: bool = true;
        for v in re {
            if cell.contains(&v) == false {
                fin = false;
                break;
            }
        }
        if fin == true {
            // println!("=>{}", theid);
            //-- check every triangle for encroachment
            let lts = dt.incident_triangles_to_vertex(*theid).unwrap();
            for t in &lts {
                // println!("t {}", t);
                if startin::geom::circumcentre_encroach_bbox(
                    &dt.get_point(t.v[0]).unwrap(),
                    &dt.get_point(t.v[1]).unwrap(),
                    &dt.get_point(t.v[2]).unwrap(),
                    &gbbox,
                ) == true
                {
                    fin = false;
                    break;
                }
            }
        }
        if fin == true {
            finpts.insert(*theid);
        }
    }
    for each in &finpts {
        let p = dt.get_point(*each).unwrap();
        io::stdout().write_all(
            &format!(
                "v {} {} {} {} {:?}\n",
                *each,
                p[0],
                p[1],
                p[2],
                dt.adjacent_vertices_to_vertex(*each).unwrap()
            )
            .as_bytes(),
        )?;
        cell.remove(each);
    }
    Ok(())
}

fn get_gx_gy(x: f64, y: f64, minx: f64, miny: f64, cellsize: usize) -> (usize, usize) {
    (
        ((x - minx) / cellsize as f64) as usize,
        ((y - miny) / cellsize as f64) as usize,
    )
}

fn get_cell_bbox(gx: usize, gy: usize, bbox: &[f64], cellsize: usize, gbbox: &mut [f64]) {
    gbbox[0] = bbox[0] + (gx * cellsize) as f64;
    gbbox[1] = bbox[1] + (gy * cellsize) as f64;
    gbbox[2] = bbox[0] + ((gx + 1) * cellsize) as f64;
    gbbox[3] = bbox[1] + ((gy + 1) * cellsize) as f64;
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
