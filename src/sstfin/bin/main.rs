use hashbrown::HashMap;
use rand::{thread_rng, Rng};
use std::env;
use std::fmt;
use std::fs::File;
use std::io::Seek;
use std::io::{self, Write};
use std::io::{BufRead, BufReader};

#[macro_use]
extern crate log; //info/debug/error

#[derive(Clone)]
pub struct Point {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point {
    pub fn new() -> Point {
        Point {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }
}

impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({}, {}, {})", self.x, self.y, self.z)
    }
}

fn main() {
    env_logger::init();

    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        error!("Input parameter not given: cellsize. Abort.");
        // eprintln!("Input parameter not given: cellsize");
        std::process::exit(1);
    }
    let f = File::open(&args[1]).expect("Unable to open file");
    let cellsize: usize = args[2].parse::<usize>().unwrap();

    //-- pass #1
    let re = pass_1(&f);
    let bbox = re.0;
    let totalpts: usize = re.1;
    info!("First pass ✅");
    info!("bbox={:?}", bbox);

    //-- sprinkler
    let mut rng = thread_rng();
    let mut sprinkled: HashMap<usize, Point> = HashMap::new();
    let nc: usize = (totalpts as f64 * 0.001) as usize; //-- TODO: what is a good value?
    for _i in 0..nc {
        sprinkled.insert(rng.gen_range(0, totalpts), Point::new());
    }
    info!("Sprinkled points: {}", sprinkled.len());

    //-- pass #2
    let mut g: Vec<Vec<usize>> = pass_2(&f, &bbox, cellsize, &mut sprinkled);
    info!("Second pass ✅");

    //-- pass #3
    let _re = pass_3(&f, &bbox, cellsize, &mut g, &sprinkled);
    info!("Third pass ✅");
}

fn pass_3(
    mut f: &File,
    bbox: &Vec<f64>,
    cellsize: usize,
    g: &mut Vec<Vec<usize>>,
    sprinkled: &HashMap<usize, Point>,
) -> io::Result<()> {
    let width: usize = ((bbox[2] - bbox[0]) / cellsize as f64).ceil() as usize;
    let height: usize = ((bbox[3] - bbox[1]) / cellsize as f64).ceil() as usize;
    let mut gpts: Vec<Vec<Vec<Point>>> = vec![vec![Vec::new(); height]; width];
    let _re = f.seek(std::io::SeekFrom::Start(0)); //-- reset to begining of the file
    let f = BufReader::new(f);
    //-- total number of points
    let mut total: usize = 0;
    for (i, _gx) in g.iter().enumerate() {
        for (j, _gy) in g[i].iter().enumerate() {
            total += g[i][j];
        }
    }
    io::stdout().write_all(b"# sstfin\n")?;
    io::stdout().write_all(&format!("n {}\n", total).as_bytes())?;
    //-- cellsize
    io::stdout().write_all(&format!("c {}\n", cellsize).as_bytes())?;
    io::stdout().write_all(&format!("d {} {}\n", width, height).as_bytes())?;
    //-- bbox
    io::stdout().write_all(&format!("b {} {}\n", bbox[0], bbox[1]).as_bytes())?;
    // io::stdout().write_all(&format!("b {:.3} {:.3}\n", bbox[0], bbox[1]).as_bytes())?;
    //-- cells that have no points
    for (i, _gx) in g.iter().enumerate() {
        for (j, _gy) in g[i].iter().enumerate() {
            if g[i][j] == 0 {
                io::stdout().write_all(&format!("x {} {}\n", i, j).as_bytes())?;
            }
        }
    }

    //-- chunker: promote these points at the top of the stream
    for (_, pt) in sprinkled {
        io::stdout().write_all(&format!("v {} {} {}\n", pt.x, pt.y, pt.z).as_bytes())?;
    }
    //-- read again the file
    for (i, l) in f.lines().enumerate() {
        let l = l.expect("Unable to read line");
        let v: Vec<f64> = l.split(' ').map(|s| s.parse().unwrap()).collect();
        let gxy: (usize, usize) = get_gx_gy(v[0], v[1], bbox[0], bbox[1], cellsize);
        // println!("{}--{}", gxy.0, gxy.1);
        g[gxy.0][gxy.1] -= 1;
        if !sprinkled.contains_key(&i) {
            gpts[gxy.0][gxy.1].push(Point {
                x: v[0],
                y: v[1],
                z: v[2],
            });
        }
        if g[gxy.0][gxy.1] == 0 {
            for pt in gpts[gxy.0][gxy.1].iter() {
                io::stdout().write_all(&format!("v {} {} {}\n", pt.x, pt.y, pt.z).as_bytes())?;
            }
            io::stdout().write_all(&format!("x {} {}\n", gxy.0, gxy.1).as_bytes())?;
            gpts[gxy.0][gxy.1].clear();
            gpts[gxy.0][gxy.1].shrink_to_fit();
        }
    }
    Ok(())
}

fn pass_2(
    mut f: &File,
    bbox: &Vec<f64>,
    cellsize: usize,
    sprinkled: &mut HashMap<usize, Point>,
) -> Vec<Vec<usize>> {
    let width: usize = ((bbox[2] - bbox[0]) / cellsize as f64).ceil() as usize;
    let height: usize = ((bbox[3] - bbox[1]) / cellsize as f64).ceil() as usize;
    info!(
        "Virtual grid is {}x{}; cellsize={}",
        width, height, cellsize
    );
    let mut g: Vec<Vec<usize>> = vec![vec![0; height]; width];
    let _re = f.seek(std::io::SeekFrom::Start(0)); //-- reset to begining of the file
    let f = BufReader::new(f);
    for (i, l) in f.lines().enumerate() {
        let l = l.expect("Unable to read line");
        let v: Vec<f64> = l.split(' ').map(|s| s.parse().unwrap()).collect();
        let gxy: (usize, usize) = get_gx_gy(v[0], v[1], bbox[0], bbox[1], cellsize);
        // println!("{}--{}", gxy.0, gxy.1);
        g[gxy.0][gxy.1] += 1;
        //-- chunking
        if sprinkled.contains_key(&i) == true {
            let pc = sprinkled.entry(i).or_insert(Point::new());
            *pc = Point {
                x: v[0],
                y: v[1],
                z: v[2],
            };
        }
    }
    g
}

fn get_gx_gy(x: f64, y: f64, minx: f64, miny: f64, cellsize: usize) -> (usize, usize) {
    (
        ((x - minx) / cellsize as f64) as usize,
        ((y - miny) / cellsize as f64) as usize,
    )
}

fn pass_1(f: &File) -> (Vec<f64>, usize) {
    let f = BufReader::new(f);
    let mut bbox: Vec<f64> = Vec::new();
    bbox.push(std::f64::MAX);
    bbox.push(std::f64::MAX);
    bbox.push(std::f64::MIN);
    bbox.push(std::f64::MIN);
    let mut n: usize = 0;
    for l in f.lines() {
        let l = l.expect("Unable to read line");
        let v: Vec<f64> = l.split(' ').map(|s| s.parse().unwrap()).collect();
        n += 1;
        //-- minxy
        for i in 0..2 {
            if v[i] < bbox[i] {
                bbox[i] = v[i];
            }
        }
        //-- maxxy
        for i in 0..2 {
            if v[i] > bbox[i + 2] {
                bbox[i + 2] = v[i];
            }
        }
    }
    (bbox, n)
}
