#[allow(dead_code)]
#[allow(unused_variables)]
extern crate startin;

use std::io::BufRead;

fn main() {
    let mut totalpts: usize = 0;
    let mut cellsize: usize = 0;
    let mut bbox: [f64; 4] = [std::f64::MAX, std::f64::MAX, std::f64::MIN, std::f64::MIN];
    let mut gpts: Vec<Vec<Vec<usize>>>;

    let mut dt = startin::Triangulation::new();

    let stdin = std::io::stdin();
    for line in stdin.lock().lines() {
        let l = line.unwrap();
        // println!("{}", l);
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
            'v' => {
                let re = parse_3_f64(&l);
                let _re = dt.insert_one_pt(re.0, re.1, re.2);
            }
            'x' => continue,
            // 'x' => println!("Cell '{}' is finalised.", l),
            _ => println!("WRONGLY FORMATTED STREAM. ABORT."),
        }
    }

    println!("\n===\nDT # points: {}", dt.number_of_vertices());
    // println!("totalpts: {}", totalpts);
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
