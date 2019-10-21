#[allow(dead_code)]
#[allow(unused_variables)]
extern crate startin;

#[macro_use]
extern crate log; //info/debug/error
use std::io::BufRead;

fn main() {
    env_logger::init();
    let mut totalpts: usize = 0;
    let mut cellsize: usize = 0;
    let mut bbox: [f64; 2] = [std::f64::MIN, std::f64::MIN];
    let mut gpts: Vec<Vec<Vec<usize>>> = Vec::new();

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
                    each.resize(re.1, Vec::new());
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
                let pos: usize;
                match re {
                    Ok(x) => pos = x,
                    Err(x) => pos = x,
                };
                // println!("999: {}", pos);
                let g = get_gx_gy(v.0, v.1, bbox[0], bbox[1], cellsize);
                gpts[g.0][g.1].push(pos);
            }
            'x' => {
                //-- finalise a cell
                let re = parse_2_usize(&l);
                finalise_cell(&mut dt, &mut gpts, (re.0, re.1));
                // println!("Finalise cell {}", l);
            }
            _ => error!("Wrongly formatted stream. Abort."),
        }
    }

    info!("DT # points: {}", dt.number_of_vertices());
}

fn finalise_cell(
    dt: &mut startin::Triangulation,
    gpts: &mut Vec<Vec<Vec<usize>>>,
    cell: (usize, usize),
) {
    // println!("Finalise_cell() {}--{}", cell.0, cell.1);
    info!(
        "cell finalised {}--{} ({} vertices)",
        cell.0,
        cell.1,
        gpts[cell.0][cell.1].len()
    );
}

fn get_gx_gy(x: f64, y: f64, minx: f64, miny: f64, cellsize: usize) -> (usize, usize) {
    (
        ((x - minx) / cellsize as f64) as usize,
        ((y - miny) / cellsize as f64) as usize,
    )
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
