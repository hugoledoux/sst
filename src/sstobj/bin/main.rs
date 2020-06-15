//-- sstobj

#[allow(dead_code)]
#[allow(unused_variables)]
#[macro_use]
extern crate log; //info/debug/error

use std::io;
use std::io::BufRead;
use std::io::Write;

use hashbrown::HashMap;

fn main() -> io::Result<()> {
    env_logger::init();

    let mut vids: HashMap<usize, usize> = HashMap::new();

    let mut vertices: Vec<String> = Vec::new();
    let mut faces: Vec<String> = Vec::new();

    info!("===== sstobj =====");

    let stdin = std::io::stdin();
    for line in stdin.lock().lines() {
        let l = line.unwrap();
        // println!("sstobj: {}", l);
        if l.is_empty() {
            continue;
        }
        let ch = l.chars().next().unwrap();
        match ch {
            '#' => continue,
            'v' => {
                let ls: Vec<&str> = l.split_whitespace().collect();
                let oid: usize = ls[1].parse::<usize>().unwrap();
                let x: f64 = ls[2].parse::<f64>().unwrap();
                let y: f64 = ls[3].parse::<f64>().unwrap();
                let z: f64 = ls[4].parse::<f64>().unwrap();
                let nid: usize = vids.len() + 1;
                vids.insert(oid, nid);
                // io::stdout().write_all(&format!("v {} {} {}\n", x, y, z).as_bytes())?;
                vertices.push(format!("v {} {} {}\n", x, y, z));
            }
            'f' => {
                let ls: Vec<&str> = l.split_whitespace().collect();
                let f1: usize = ls[1].parse::<usize>().unwrap();
                let f2: usize = ls[2].parse::<usize>().unwrap();
                let f3: usize = ls[3].parse::<usize>().unwrap();
                // vids.insert(oid, nid);
                // io::stdout().write_all(
                //     &format!("f {} {} {}\n", vids[&f1], vids[&f2], vids[&f3]).as_bytes(),
                // )?;
                faces.push(format!("f {} {} {}\n", vids[&f1], vids[&f2], vids[&f3]));
            }
            'x' => {
                let ls: Vec<&str> = l.split_whitespace().collect();
                let _v: usize = ls[1].parse::<usize>().unwrap();
                // info!("vertex #{} finalised", v);
                // vids.remove(&v);
            }

            _ => {
                error!("Wrongly formatted stream. Abort.");
                std::process::exit(1);
            }
        }
    }
    info!("Finished reading the stream");

    for v in vertices {
        io::stdout().write_all(v.as_bytes())?;
    }
    for f in faces {
        io::stdout().write_all(f.as_bytes())?;
    }
    Ok(())
}
