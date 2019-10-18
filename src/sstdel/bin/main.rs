#[allow(dead_code)]
#[allow(unused_variables)]
extern crate startin;

use std::io::BufRead;

fn main() {
    let mut totalpts: usize = 0;
    let mut cellsize: usize = 0;
    let mut bbox: [f64; 4] = [std::f64::MAX, std::f64::MAX, std::f64::MIN, std::f64::MIN];
    let mut gpts: Vec<Vec<Vec<usize>>>;

    let stdin = std::io::stdin();
    for line in stdin.lock().lines() {
        let l = line.unwrap();
        if l.is_empty() {
            continue;
        }
        let ch = l.chars().next().unwrap();
        match ch {
            '#' => continue,
            'n' => {
                totalpts = l
                    .split_whitespace()
                    .last()
                    .unwrap()
                    .parse::<usize>()
                    .unwrap()
            }
            'c' => {
                cellsize = l
                    .split_whitespace()
                    .last()
                    .unwrap()
                    .parse::<usize>()
                    .unwrap()
            }
            'd' => {
                let re = parse_2_usize(&l);
                gpts = vec![vec![Vec::new(); re.0]; re.1];
            }
            'b' => {
                let re = parse_2_f64(&l);
                bbox[0] = re.0;
                bbox[1] = re.1;
            }
            'v' => continue,
            'x' => continue,
            _ => println!("WRONGLY FORMATTED STREAM. ABORT."),
        }
    }

    let mut dt = startin::Triangulation::new();
    let _re = dt.insert_one_pt(1.0, 2.0, 3.0);
    println!("# points: {}", dt.number_of_vertices());
    println!("totalpts: {}", totalpts);
}

// fn parse_grid_dimensions(l: &String, gpts: &mut Vec<Vec<Vec<usize>>>) {
//     let mut it = l.split_whitespace();
//     let wit = it.next();
//     let hit = it.next();
//     let width: usize = wit.unwrap().parse::<usize>().unwrap();
//     it.next();
//     let height: usize = hit.unwrap().parse::<usize>().unwrap();
//     gpts = vec![vec![Vec::new(); height]; width];
// }

fn parse_2_usize(l: &String) -> (usize, usize) {
    let mut it = l.split_whitespace();
    let a: usize = it.next().unwrap().parse::<usize>().unwrap();
    let b: usize = it.next().unwrap().parse::<usize>().unwrap();
    (a, b)
}

fn parse_2_f64(l: &String) -> (f64, f64) {
    let mut it = l.split_whitespace();
    let a: f64 = it.next().unwrap().parse::<f64>().unwrap();
    let b: f64 = it.next().unwrap().parse::<f64>().unwrap();
    (a, b)
}

fn doiets() {
    println!("iets.");
}

// let l = l.expect("Unable to read line");
// let v: Vec<f64> = l.split(' ').map(|s| s.parse().unwrap()).collect();
// let gxy: (usize, usize) = get_gx_gy(v[0], v[1], bbox[0], bbox[1], cellsize);
