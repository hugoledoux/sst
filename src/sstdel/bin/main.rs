extern crate startin;
use std::io::BufRead;

fn main() {
    let stdin = std::io::stdin();
    for line in stdin.lock().lines() {
        // println!("{}", line.unwrap());
        let ch = line.unwrap().chars().next().unwrap();
        match ch {
            // '#' => ,
            'c' => doiets(),
            _ => println!("WRONGLY FORMATTED STREAM. ABORT."),
        }
        // if ch == '#' {
        //     println!("ok");
        // }
    }

    let mut dt = startin::Triangulation::new();
    let _re = dt.insert_one_pt(1.0, 2.0, 3.0);

    println!("# points: {}", dt.number_of_vertices());
}

fn doiets() {
    println!("iets.");
}
