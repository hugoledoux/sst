//-- sstobj

#[allow(dead_code)]
#[allow(unused_variables)]
#[macro_use]
extern crate log; //info/debug/error

use std::io;
use std::io::BufRead;

fn main() -> io::Result<()> {
    env_logger::init();

    info!("Init DT");

    let stdin = std::io::stdin();
    for line in stdin.lock().lines() {
        // let l = l.expect("Unable to read line");
        let l = line.unwrap();
        println!("sstobj: {}", l);
        if l.is_empty() {
            continue;
        }
        // let ch = l.chars().next().unwrap();
        // match ch {
        //     '#' => continue,
        //     'v' => {
        //         _totalpts = l
        //             .split_whitespace()
        //             .last()
        //             .unwrap()
        //             .parse::<usize>()
        //             .unwrap()
        //     }
        //     _ => {
        //         error!("Wrongly formatted stream. Abort.");
        //         std::process::exit(1);
        //     }
        // }
    }
    info!("Finished reading the stream");

    Ok(())
}
