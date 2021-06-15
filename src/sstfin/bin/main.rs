//-- sstfin

use hashbrown::HashMap;
use rand::{thread_rng, Rng};

use std::fmt;
use std::fs::File;

use std::io::{self, Write};
use std::io::{BufRead, BufReader};
use std::path::Path;

use clap::App;
use num_format::{Locale, ToFormattedString};

extern crate las;

use las::Read;

#[macro_use]
extern crate log; //info/debug/error

#[derive(Clone)]
pub struct Point {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(PartialEq, Eq)]
enum InputType {
    LAS,
    XYZ,
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
    let matches = App::new("sstfin")
        .version("0.1")
        .about("streaming startin -- finalisation")
        .arg("<INPUT>             'The input LAS/LAZ/XYZ file(s) to use.'")
        .arg("<RESOLUTION>        'The cell resolution (integer, eg '5' for 5mX5m)'")
        .arg("-g...               'Use only ground in LAS files'")
        .arg("--sprinkle...       'Value to use (0.001 is default; totalpts * 0.001)'")
        .get_matches();

    env_logger::init();

    let mut inputformat: InputType = InputType::LAS;

    let mut paths: Vec<String> = Vec::new();
    let t1 = matches.value_of("INPUT").unwrap();
    let path = Path::new(&t1);
    if path.extension().unwrap() == "las"
        || path.extension().unwrap() == "laz"
        || path.extension().unwrap() == "LAS"
        || path.extension().unwrap() == "LAZ"
    {
        paths.push(path.to_str().unwrap().to_string());
    } else if path.extension().unwrap() == "xyz" {
        inputformat = InputType::XYZ;
        paths.push(path.to_str().unwrap().to_string());
    } else if path.extension().unwrap() == "files" {
        let f = File::open(path).expect("Unable to open file");
        let f = BufReader::new(f);
        for l in f.lines() {
            let l = l.expect("Unable to read line");
            paths.push(l);
        }
    } else {
        eprintln!("input format not accepted");
        std::process::exit(1);
    }

    info!("Number of input files: {}", paths.len());
    for (i, each) in paths.iter().enumerate() {
        info!("\t{}. {}", i + 1, each);
    }

    let cellsize: usize = matches.value_of("RESOLUTION").unwrap().parse().unwrap();

    let mut sprinkle_param: f64 = 0.001;
    if matches.occurrences_of("sprinkle") > 0 {
        sprinkle_param = matches.value_of("sprinkle").unwrap().parse().unwrap();
    }

    //-- pass #1
    let re = pass_1(&paths, &inputformat);
    let bbox = re.0;
    let totalpts: usize = re.1;
    info!("count pts={}", totalpts.to_formatted_string(&Locale::en));
    info!("bbox={:?}", bbox);
    info!("First pass 🏁");

    //-- sprinkler
    let mut rng = thread_rng();
    let mut sprinkled: HashMap<usize, Point> = HashMap::new();
    let nc: usize = (totalpts as f64 * sprinkle_param) as usize; //-- TODO: what is a good value?
    for _i in 0..nc {
        sprinkled.insert(rng.gen_range(0..totalpts), Point::new());
    }
    info!("Sprinkled points: {}", sprinkled.len());

    //-- pass #2
    info!("Second pass ➡️");
    let mut g: Vec<Vec<usize>> = pass_2(&paths, &bbox, cellsize, &mut sprinkled, &inputformat);
    info!("Second pass ✅");

    // //-- pass #3
    info!("Third pass ➡️");
    let _re = pass_3(&paths, &bbox, cellsize, &mut g, &sprinkled, &inputformat);
    info!("Third pass ✅");
}

fn pass_1(paths: &Vec<String>, it: &InputType) -> (Vec<f64>, usize) {
    if *it == InputType::LAS {
        pass_1_las(&paths)
    } else {
        pass_1_xyz(&paths)
    }
}

fn pass_1_las(paths: &Vec<String>) -> (Vec<f64>, usize) {
    let mut n: usize = 0;
    let mut bbox = vec![std::f64::MAX, std::f64::MAX, std::f64::MIN, std::f64::MIN];
    for path in paths {
        let reader = las::Reader::from_path(path).expect("Unable to open reader");
        n += reader.header().number_of_points() as usize;
        let b = reader.header().bounds();
        if b.min.x < bbox[0] {
            bbox[0] = b.min.x;
        }
        if b.min.y < bbox[1] {
            bbox[1] = b.min.y;
        }
        if b.max.x > bbox[2] {
            bbox[2] = b.max.x;
        }
        if b.max.y > bbox[3] {
            bbox[3] = b.max.y;
        }
    }
    (bbox, n)
}

fn pass_1_xyz(paths: &Vec<String>) -> (Vec<f64>, usize) {
    let mut n: usize = 0;
    let mut bbox = vec![std::f64::MAX, std::f64::MAX, std::f64::MIN, std::f64::MIN];
    for path in paths {
        let f = File::open(path).expect("Unable to open file");
        let f = BufReader::new(f);
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
    }
    (bbox, n)
}

fn pass_2(
    paths: &Vec<String>,
    bbox: &Vec<f64>,
    cellsize: usize,
    mut sprinkled: &mut HashMap<usize, Point>,
    it: &InputType,
) -> Vec<Vec<usize>> {
    if *it == InputType::LAS {
        pass_2_las(&paths, &bbox, cellsize, &mut sprinkled)
    } else {
        pass_2_xyz(&paths, &bbox, cellsize, &mut sprinkled)
    }
}

fn pass_2_las(
    paths: &Vec<String>,
    bbox: &Vec<f64>,
    cellsize: usize,
    sprinkled: &mut HashMap<usize, Point>,
) -> Vec<Vec<usize>> {
    let width = (bbox[2] - bbox[0]) / cellsize as f64;
    let height = (bbox[3] - bbox[1]) / cellsize as f64;
    //-- make it a square to have a quadtree
    let mut tmp = width;
    if height > width {
        tmp = height;
    }
    //-- needs to be a power^2 so that is quadtree
    let cellno: u32 = (2_u32).pow(tmp.log(2.0).ceil() as u32);
    info!(
        "Virtual grid is {}x{}; cellsize={}",
        cellno, cellno, cellsize
    );
    let mut g: Vec<Vec<usize>> = vec![vec![0; cellno as usize]; cellno as usize];

    let mut count: usize = 0;
    for path in paths {
        let mut reader = las::Reader::from_path(path).expect("Unable to open reader");
        for each in reader.points() {
            let p = each.unwrap();
            let gxy: (usize, usize) = get_gx_gy(p.x, p.y, bbox[0], bbox[1], cellsize);
            g[gxy.0][gxy.1] += 1;
            //-- chunking
            if sprinkled.contains_key(&count) == true {
                let pc = sprinkled.entry(count).or_insert(Point::new());
                *pc = Point {
                    x: p.x,
                    y: p.y,
                    z: p.z,
                };
            }
            count += 1;
        }
    }
    g
}

fn pass_2_xyz(
    paths: &Vec<String>,
    bbox: &Vec<f64>,
    cellsize: usize,
    sprinkled: &mut HashMap<usize, Point>,
) -> Vec<Vec<usize>> {
    let width = (bbox[2] - bbox[0]) / cellsize as f64;
    let height = (bbox[3] - bbox[1]) / cellsize as f64;
    //-- make it a square to have a quadtree
    let mut tmp = width;
    if height > width {
        tmp = height;
    }
    //-- needs to be a power^2 so that is quadtree
    let cellno: u32 = (2_u32).pow(tmp.log(2.0).ceil() as u32);
    info!(
        "Virtual grid is {}x{}; cellsize={}",
        cellno, cellno, cellsize
    );
    let mut g: Vec<Vec<usize>> = vec![vec![0; cellno as usize]; cellno as usize];
    let mut count: usize = 0;
    for path in paths {
        let f = File::open(path).expect("Unable to open file");
        let f = BufReader::new(f);
        for l in f.lines() {
            let l = l.expect("Unable to read line");
            let v: Vec<f64> = l.split(' ').map(|s| s.parse().unwrap()).collect();
            let gxy: (usize, usize) = get_gx_gy(v[0], v[1], bbox[0], bbox[1], cellsize);
            g[gxy.0][gxy.1] += 1;
            //-- chunking
            if sprinkled.contains_key(&count) == true {
                let pc = sprinkled.entry(count).or_insert(Point::new());
                *pc = Point {
                    x: v[0],
                    y: v[1],
                    z: v[2],
                };
            }
            count += 1;
        }
    }
    g
}

fn pass_3(
    paths: &Vec<String>,
    bbox: &Vec<f64>,
    cellsize: usize,
    mut g: &mut Vec<Vec<usize>>,
    sprinkled: &HashMap<usize, Point>,
    it: &InputType,
) -> io::Result<()> {
    if *it == InputType::LAS {
        pass_3_las(&paths, &bbox, cellsize, &mut g, &sprinkled)
    } else {
        pass_3_xyz(&paths, &bbox, cellsize, &mut g, &sprinkled)
    }
}

fn pass_3_las(
    paths: &Vec<String>,
    bbox: &Vec<f64>,
    cellsize: usize,
    g: &mut Vec<Vec<usize>>,
    sprinkled: &HashMap<usize, Point>,
) -> io::Result<()> {
    let cellno: usize = g[0].len();
    let mut gpts: Vec<Vec<Vec<Point>>> = vec![vec![Vec::new(); cellno]; cellno];

    //-- total number of points
    let mut total: usize = 0;
    for (i, _gx) in g.iter().enumerate() {
        for (j, _gy) in g[i].iter().enumerate() {
            total += g[i][j];
        }
    }
    io::stdout().write_all(b"# sstfin\n")?;
    io::stdout().write_all(&format!("n {}\n", total).as_bytes())?;
    //-- number of cells cXc
    io::stdout().write_all(&format!("c {}\n", cellno).as_bytes())?;
    //-- cellsize
    io::stdout().write_all(&format!("s {}\n", cellsize).as_bytes())?;
    //-- bbox
    io::stdout().write_all(
        &format!(
            "b {:.3} {:.3} {:.3} {:.3}\n",
            bbox[0], bbox[1], bbox[2], bbox[3]
        )
        .as_bytes(),
    )?;

    //-- cells that have no points
    for (i, _gx) in g.iter().enumerate() {
        for (j, _gy) in g[i].iter().enumerate() {
            if g[i][j] == 0 {
                io::stdout().write_all(&format!("x {} {}\n", i, j).as_bytes())?;
            }
        }
    }

    //-- find the transform.scale
    let reader = las::Reader::from_path(&paths[0]).expect("Unable to open reader");
    let tmp = reader.header().transforms().x.scale;
    let impdigits = tmp.log10().abs().round() as usize;
    info!("impdigits {}", impdigits);

    //-- chunker: promote these points at the top of the stream
    for (_, pt) in sprinkled {
        io::stdout().write_all(
            &format!("v {0:.3$} {1:.3$} {2:.3$}\n", pt.x, pt.y, pt.z, impdigits).as_bytes(),
        )?;
    }

    //-- read again the files
    let mut count: usize = 0;
    for path in paths {
        let mut reader = las::Reader::from_path(path).expect("Unable to open reader");
        for each in reader.points() {
            let p = each.unwrap();
            let gxy: (usize, usize) = get_gx_gy(p.x, p.y, bbox[0], bbox[1], cellsize);
            g[gxy.0][gxy.1] -= 1;
            if !sprinkled.contains_key(&count) {
                gpts[gxy.0][gxy.1].push(Point {
                    x: p.x,
                    y: p.y,
                    z: p.z,
                });
            }
            if g[gxy.0][gxy.1] == 0 {
                for pt in gpts[gxy.0][gxy.1].iter() {
                    io::stdout().write_all(
                        &format!("v {0:.3$} {1:.3$} {2:.3$}\n", pt.x, pt.y, pt.z, impdigits)
                            .as_bytes(),
                    )?;
                }
                io::stdout().write_all(&format!("x {} {}\n", gxy.0, gxy.1).as_bytes())?;
                gpts[gxy.0][gxy.1].clear();
                gpts[gxy.0][gxy.1].shrink_to_fit();
            }
            count += 1;
        }
    }
    Ok(())
}

fn pass_3_xyz(
    paths: &Vec<String>,
    bbox: &Vec<f64>,
    cellsize: usize,
    g: &mut Vec<Vec<usize>>,
    sprinkled: &HashMap<usize, Point>,
) -> io::Result<()> {
    // let width: usize = ((bbox[2] - bbox[0]) / cellsize as f64).ceil() as usize;
    // let height: usize = ((bbox[3] - bbox[1]) / cellsize as f64).ceil() as usize;
    let cellno: usize = g[0].len();
    let mut gpts: Vec<Vec<Vec<Point>>> = vec![vec![Vec::new(); cellno]; cellno];

    //-- total number of points
    let mut total: usize = 0;
    for (i, _gx) in g.iter().enumerate() {
        for (j, _gy) in g[i].iter().enumerate() {
            total += g[i][j];
        }
    }
    io::stdout().write_all(b"# sstfin\n")?;
    io::stdout().write_all(&format!("n {}\n", total).as_bytes())?;
    //-- number of cells cXc
    io::stdout().write_all(&format!("c {}\n", cellno).as_bytes())?;
    //-- cellsize
    io::stdout().write_all(&format!("s {}\n", cellsize).as_bytes())?;
    //-- bbox
    io::stdout().write_all(
        &format!(
            "b {:.3} {:.3} {:.3} {:.3}\n",
            bbox[0], bbox[1], bbox[2], bbox[3]
        )
        .as_bytes(),
    )?;

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

    //-- read again the files
    //-- read again the files
    let mut count: usize = 0;
    for path in paths {
        let f = File::open(path).expect("Unable to open file");
        let f = BufReader::new(f);
        for l in f.lines() {
            let l = l.expect("Unable to read line");
            let v: Vec<f64> = l.split(' ').map(|s| s.parse().unwrap()).collect();
            let gxy: (usize, usize) = get_gx_gy(v[0], v[1], bbox[0], bbox[1], cellsize);
            g[gxy.0][gxy.1] -= 1;
            if !sprinkled.contains_key(&count) {
                gpts[gxy.0][gxy.1].push(Point {
                    x: v[0],
                    y: v[1],
                    z: v[2],
                });
            }
            if g[gxy.0][gxy.1] == 0 {
                for pt in gpts[gxy.0][gxy.1].iter() {
                    io::stdout()
                        .write_all(&format!("v {} {} {}\n", pt.x, pt.y, pt.z).as_bytes())?;
                }
                io::stdout().write_all(&format!("x {} {}\n", gxy.0, gxy.1).as_bytes())?;
                gpts[gxy.0][gxy.1].clear();
                gpts[gxy.0][gxy.1].shrink_to_fit();
            }
            count += 1;
        }
    }
    Ok(())
}

fn get_gx_gy(x: f64, y: f64, minx: f64, miny: f64, cellsize: usize) -> (usize, usize) {
    (
        ((x - minx) / cellsize as f64) as usize,
        ((y - miny) / cellsize as f64) as usize,
    )
}
