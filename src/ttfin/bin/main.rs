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
    /// Use only ground points (LAS classification == 2)
    #[arg(short, long)]
    ground: bool,
    /// Value to use (totalpts * 0.001)
    #[arg(short, long, default_value_t = 0.001)]
    sprinkle: f64,
    /// Write the ASCII `.spa` (default is the binary format)
    #[arg(long)]
    ascii: bool,
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

    if cli.ground {
        info!("Using only ground points (classification == 2)");
    }

    //-- pass #1
    let re = pass_1(&paths);
    let bbox = re.0;
    let totalpts: usize = re.1;
    if cli.ground == false {
        info!("count pts={}", totalpts.to_formatted_string(&Locale::en));
    }
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
    let mut g: Vec<Vec<usize>> = pass_2(&paths, &bbox, cellsize, &mut sprinkled, cli.ground);
    if cli.ground {
        //-- the header count includes every class; report how many ground points
        //-- are actually kept (this is the `n` that will be written out)
        let kept: usize = g.iter().flatten().sum();
        info!("count pts={} (ground only)", kept.to_formatted_string(&Locale::en));
    }
    info!("Second pass 🏁");

    //-- the CRS of the input (E)VLRs, carried through the stream as a `p`
    //-- record so that downstream tools (eg ttdem) can georeference it
    let crs = extract_crs(&paths);
    if crs.0 != 0 {
        info!("CRS: EPSG:{} (vertical EPSG:{})", crs.0, crs.1);
    } else {
        info!("no CRS found in the input");
    }

    // //-- pass #3
    info!("Third pass ➡️");
    let _re = pass_3(&paths, &bbox, cellsize, &mut g, &sprinkled, cli.ground, cli.ascii == false, crs);
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

/// Extract (horizontal, vertical) EPSG codes from the CRS (E)VLRs of the
/// first input file, to emit as a `p` record in the stream.
///
/// GeoTIFF keys VLRs give the codes directly; a WKT VLR (or LAS 1.4 WKT
/// EVLR) is the fallback. In a WKT the CRS's own AUTHORITY is the last
/// EPSG code inside its element, so the horizontal code comes from the
/// PROJCS/PROJCRS slice and the vertical one from the VERT_CS/VERTCRS slice
/// (the nested GEOGCS/VERT_DATUM authorities that appear earlier in those
/// slices are ignored this way).
fn extract_crs(paths: &Vec<String>) -> (u32, u32) {
    let reader = match las::Reader::from_path(&paths[0]) {
        Ok(r) => r,
        Err(e) => {
            info!("could not reopen {} for the CRS: {}", paths[0], e);
            return (0, 0);
        }
    };
    let header = reader.header();

    if let Ok(Some(gk)) = header.get_geotiff_crs() {
        let h = gk.get_projected_crs_geo_key_value().unwrap_or(0);
        let v = gk.get_vertical_crs_geo_key_value().unwrap_or(0);
        if h != 0 {
            return (h as u32, v as u32);
        }
    }

    let Some(wkt) = header
        .get_wkt_crs_bytes()
        .map(|b| String::from_utf8_lossy(b).into_owned())
    else {
        info!("no CRS in the input; the output GeoTIFF will have no CRS");
        return (0, 0);
    };
    let h = wkt_element_slice(&wkt, "PROJCS")
        .or_else(|| wkt_element_slice(&wkt, "PROJCRS"))
        .as_deref()
        .and_then(|s| epsg_codes_in_wkt(s).pop())
        .or_else(|| epsg_codes_in_wkt(&wkt).pop()) //-- plain GEOGCS
        .unwrap_or(0);
    let v = wkt_element_slice(&wkt, "VERT_CS")
        .or_else(|| wkt_element_slice(&wkt, "VERTCRS"))
        .as_deref()
        .and_then(|s| epsg_codes_in_wkt(s).pop())
        .unwrap_or(0);
    (h, v)
}

/// The slice of a WKT from the opening bracket of `keyword` to its matching
/// close bracket, so nested elements can be scanned in isolation.
fn wkt_element_slice(wkt: &str, keyword: &str) -> Option<String> {
    let start = wkt.find(keyword)? + keyword.len();
    let open = wkt[start..].find(['[', '('])? + start;
    let mut depth = 0usize;
    for (k, b) in wkt.bytes().enumerate().skip(open) {
        match b {
            b'[' | b'(' => depth += 1,
            b']' | b')' => {
                depth -= 1;
                if depth == 0 {
                    return Some(wkt[open..k + 1].to_string());
                }
            }
            _ => {}
        }
    }
    None
}

/// Every `AUTHORITY["EPSG", "<code>"]` (WKT1) or `ID["EPSG", <code>]` (WKT2)
/// code in a WKT string, in order of appearance.
fn epsg_codes_in_wkt(wkt: &str) -> Vec<u32> {
    let mut codes = Vec::new();
    let mut rest = wkt;
    while let Some(i) = rest.find("EPSG") {
        rest = &rest[i + 4..];
        let digits: String = rest
            .trim_start_matches(|c: char| !c.is_ascii_digit())
            .chars()
            .take_while(|c| c.is_ascii_digit())
            .collect();
        if let Ok(c) = digits.parse::<u32>() {
            codes.push(c);
        }
    }
    codes
}

fn pass_2(
    paths: &Vec<String>,
    bbox: &Vec<f64>,
    cellsize: usize,
    sprinkled: &mut HashMap<usize, Point>,
    ground: bool,
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
                //-- chunking: fill the sprinkled coordinates *before* the
                //-- ground filter below; with -g, a sprinkled index on a
                //-- non-ground point must be dropped here, otherwise its
                //-- (0,0,0) default would be emitted at the top of the stream
                if ground && p.classification != las::point::Classification::Ground {
                    sprinkled.remove(&count);
                } else if let Some(pc) = sprinkled.get_mut(&count) {
                    *pc = Point {
                        x: p.x,
                        y: p.y,
                        z: p.z,
                    };
                }
                //-- keep `count` running over *all* points so that the
                //-- sprinkled indices stay aligned with the header count
                if ground && p.classification != las::point::Classification::Ground {
                    count += 1;
                    continue;
                }
                let gxy: (usize, usize) = get_gx_gy(p.x, p.y, bbox[0], bbox[1], cellsize);
                g[gxy.0][gxy.1] += 1;
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
    ground: bool,
    binary: bool,
    crs: (u32, u32),
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
    if binary == true {
        sst::stream::spa_write_magic(&mut w)?;
        sst::stream::spa_write(&mut w, &sst::stream::SpaRecord::N(total as u64))?;
        sst::stream::spa_write(&mut w, &sst::stream::SpaRecord::C(cellno as u32))?;
        sst::stream::spa_write(&mut w, &sst::stream::SpaRecord::S(cellsize as u64))?;
        sst::stream::spa_write(
            &mut w,
            &sst::stream::SpaRecord::B([bbox[0], bbox[1], bbox[2], bbox[3]]),
        )?;
        sst::stream::spa_write(&mut w, &sst::stream::SpaRecord::P(crs.0, crs.1))?;
    } else {
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
        //-- CRS: horizontal and vertical EPSG (0 = unknown)
        write!(w, "p {} {}\n", crs.0, crs.1)?;
    }

    //-- cells that have no points
    for (i, _gx) in g.iter().enumerate() {
        for (j, _gy) in g[i].iter().enumerate() {
            if g[i][j] == 0 {
                if binary == true {
                    sst::stream::spa_write(&mut w, &sst::stream::SpaRecord::X(i as u32, j as u32))?;
                } else {
                    write!(w, "x {} {}\n", i, j)?;
                }
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
        if binary == true {
            sst::stream::spa_write(&mut w, &sst::stream::SpaRecord::W(pt.x, pt.y, pt.z))?;
        } else {
            write!(
                w,
                "w {0:.3$} {1:.3$} {2:.3$}\n",
                pt.x, pt.y, pt.z, impdigits
            )?;
        }
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
                //-- same filter, and same `count` treatment, as pass_2
                if ground && p.classification != las::point::Classification::Ground {
                    count += 1;
                    continue;
                }
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
                        if binary == true {
                            sst::stream::spa_write(
                                &mut w,
                                &sst::stream::SpaRecord::V(pt.x, pt.y, pt.z),
                            )?;
                        } else {
                            write!(
                                w,
                                "v {0:.3$} {1:.3$} {2:.3$}\n",
                                pt.x, pt.y, pt.z, impdigits
                            )?;
                        }
                    }
                    if binary == true {
                        sst::stream::spa_write(
                            &mut w,
                            &sst::stream::SpaRecord::X(gxy.0 as u32, gxy.1 as u32),
                        )?;
                    } else {
                        write!(w, "x {} {}\n", gxy.0, gxy.1)?;
                    }
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
