//-- ttfin

use hashbrown::HashMap;
use rand::RngExt;

use std::fmt;
use std::fs::File;

use std::io::{self, Write};
use std::io::{BufRead, BufReader};
use std::path::Path;

use num_format::{Locale, ToFormattedString};

extern crate las;

/// number of points decoded per chunk when streaming a LAS/LAZ file
const CHUNK_PTS: u64 = 1_000_000;

#[macro_use]
extern crate log; //info/debug/error

use clap::Parser;
#[derive(Parser)]
#[command(name = "ttfin")]
#[command(about = "tintin -- finalisation [ttfin]")]
#[command(author, version)]
struct Cli {
    /// The input LAS/LAZ/XYZ file(s) to use
    input: String,
    /// The cell resolution (integer, eg '5' for 5mX5m)
    resolution: usize,
    /// Use only ground in LAS files
    #[arg(short, long)]
    ground: bool,
    /// Value to use (totalpts * 0.001)
    #[arg(short, long, default_value_t = 0.001)]
    sprinkle: f64,
    #[clap(flatten)]
    verbose: clap_verbosity_flag::Verbosity,
}

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
    let cli = Cli::parse();
    env_logger::Builder::new()
        .filter_level(cli.verbose.log_level_filter())
        .init();

    let mut paths: Vec<String> = Vec::new();
    let path = Path::new(&cli.input);
    if path.extension().unwrap() == "las"
        || path.extension().unwrap() == "laz"
        || path.extension().unwrap() == "LAS"
        || path.extension().unwrap() == "LAZ"
    {
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

    let cellsize: usize = cli.resolution;

    let sprinkle_param = cli.sprinkle;

    //-- pass #1
    let re = pass_1(&paths);
    let bbox = re.0;
    let totalpts: usize = re.1;
    info!("count pts={}", totalpts.to_formatted_string(&Locale::en));
    info!("bbox={:?}", bbox);
    info!("First pass 🏁");

    //-- sprinkler
    let mut rng = rand::rng();
    let mut sprinkled: HashMap<usize, Point> = HashMap::new();
    let nc: usize = (totalpts as f64 * sprinkle_param) as usize; //-- TODO: what is a good value?
    for _i in 0..nc {
        sprinkled.insert(rng.random_range(0..totalpts), Point::new());
    }
    info!("Sprinkled points: {}", sprinkled.len());

    //-- pass #2
    info!("Second pass ➡️");
    let mut g: Vec<Vec<usize>> = pass_2(&paths, &bbox, cellsize, &mut sprinkled);
    info!("Second pass 🏁");

    // //-- pass #3
    info!("Third pass ➡️");
    let _re = pass_3(&paths, &bbox, cellsize, &mut g, &sprinkled);
    info!("Third pass 🏁");

    info!("✅");
}

fn pass_1(paths: &Vec<String>) -> (Vec<f64>, usize) {
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

fn pass_2(
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
        let mut points = las::PointDataBuilder::new()
            .for_header(reader.header())
            .build();
        loop {
            let n = reader
                .fill_points(CHUNK_PTS, &mut points)
                .expect("Unable to read points");
            if n == 0 {
                break;
            }
            for each in points.points() {
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
    }
    g
}

fn pass_3(
    paths: &Vec<String>,
    bbox: &Vec<f64>,
    cellsize: usize,
    g: &mut Vec<Vec<usize>>,
    sprinkled: &HashMap<usize, Point>,
) -> io::Result<()> {
    let cellno: usize = g[0].len();
    let mut gpts: Vec<Vec<Vec<Point>>> = vec![vec![Vec::new(); cellno]; cellno];

    //-- buffered output: one big write per chunk instead of one per line
    let stdout = io::stdout();
    let mut w = std::io::BufWriter::with_capacity(1 << 20, stdout.lock());

    //-- total number of points
    let mut total: usize = 0;
    for (i, _gx) in g.iter().enumerate() {
        for (j, _gy) in g[i].iter().enumerate() {
            total += g[i][j];
        }
    }
    w.write_all(b"# ttfin\n")?;
    write!(w, "n {}\n", total)?;
    //-- number of cells cXc
    write!(w, "c {}\n", cellno)?;
    //-- cellsize
    write!(w, "s {}\n", cellsize)?;
    //-- bbox
    write!(
        w,
        "b {:.3} {:.3} {:.3} {:.3}\n",
        bbox[0], bbox[1], bbox[2], bbox[3]
    )?;

    //-- cells that have no points
    for (i, _gx) in g.iter().enumerate() {
        for (j, _gy) in g[i].iter().enumerate() {
            if g[i][j] == 0 {
                write!(w, "x {} {}\n", i, j)?;
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
        write!(
            w,
            "w {0:.3$} {1:.3$} {2:.3$}\n",
            pt.x, pt.y, pt.z, impdigits
        )?;
    }
    w.flush()?;

    //-- read again the files
    let mut count: usize = 0;
    for path in paths {
        let mut reader = las::Reader::from_path(path).expect("Unable to open reader");
        let mut points = las::PointDataBuilder::new()
            .for_header(reader.header())
            .build();
        loop {
            let n = reader
                .fill_points(CHUNK_PTS, &mut points)
                .expect("Unable to read points");
            if n == 0 {
                break;
            }
            for each in points.points() {
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
                        write!(
                            w,
                            "v {0:.3$} {1:.3$} {2:.3$}\n",
                            pt.x, pt.y, pt.z, impdigits
                        )?;
                    }
                    write!(w, "x {} {}\n", gxy.0, gxy.1)?;
                    gpts[gxy.0][gxy.1].clear();
                    gpts[gxy.0][gxy.1].shrink_to_fit();
                }
                count += 1;
            }
        }
    }
    w.flush()?;
    Ok(())
}

fn get_gx_gy(x: f64, y: f64, minx: f64, miny: f64, cellsize: usize) -> (usize, usize) {
    (
        ((x - minx) / cellsize as f64) as usize,
        ((y - miny) / cellsize as f64) as usize,
    )
}
