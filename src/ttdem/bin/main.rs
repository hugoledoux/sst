//-- ttdem
//--
//-- Rasterises the `.sma` TIN stream of ttdt onto a regular grid by linear
//-- interpolation (Isenburg et al. 2006, §3.1), and writes the result as a
//-- Cloud Optimised GeoTIFF. Reads the stream from stdin; the COG goes to a
//-- file (TIFF needs a seekable output, so this is the one tintin operator
//-- that does not stream to stdout).
//--
//-- Grid: rows run north-to-south (row 0 = north), cell *centres* are
//-- sampled, and the bbox is the outer edge of the pixel areas. The bbox
//-- snaps outward onto a global grid aligned at --origin, so runs over
//-- neighbouring datasets line up.

mod grid;
mod raster;
mod store;
mod writer;

#[macro_use]
extern crate log; //info/debug/error

use std::io::{self, BufRead, BufReader};
use std::path::PathBuf;
use std::time::Instant;

use clap::Parser;
use geotiff_writer::{Compression, Predictor};
use hashbrown::HashMap;

use grid::Grid;
use raster::RasterStats;
use sst::stream::SmaRecord;
use store::TileStore;

/// "no data" value of the output raster; deliberately finite so that the
/// overview averaging can skip it (NaN == NaN is false).
const NODATA: f32 = -9999.0;

#[derive(Parser)]
#[command(name = "ttdem")]
#[command(about = "tintin -- rasterise the .sma stream of ttdt to a COG [ttdem]")]
#[command(author, version)]
struct Cli {
    /// Output GeoTIFF (Cloud Optimised); a seekable file, not a pipe
    output: PathBuf,
    /// Output post spacing in CRS units (metres for a projected CRS)
    #[arg(short, long)]
    res: f64,
    /// Raster bbox as outer pixel edges: minx miny maxx maxy
    /// (default: the stream bbox, snapped outward)
    #[arg(long, num_args = 4)]
    bbox: Option<Vec<f64>>,
    /// Origin of the global grid that the bbox snaps onto
    #[arg(long, num_args = 2, default_values_t = vec![0.0, 0.0])]
    origin: Vec<f64>,
    /// Horizontal CRS EPSG code ("28992" or "EPSG:28992"); overrides the
    /// stream's `p` record
    #[arg(long)]
    crs: Option<String>,
    /// Vertical CRS EPSG code (5709 for NAP); overrides the stream
    #[arg(long)]
    vertical_crs: Option<String>,
    /// Reject triangles with a 2D edge longer than this many CRS units
    /// (0 = keep all; gaps then get interpolated across)
    #[arg(long, default_value_t = 0.0)]
    max_edge: f64,
    /// COG tile size in pixels (multiple of 16)
    #[arg(long, default_value_t = 512)]
    tile_size: u32,
    /// Tile compression: deflate | lzw | none
    #[arg(long, default_value = "deflate")]
    compress: String,
    /// Overviews: auto | none
    #[arg(long, default_value = "auto")]
    overviews: String,
    /// Memory budget for the raster cache, eg 2G (the raster spills to disk
    /// beyond it; the vertex table of the stream is on top of this)
    #[arg(long, default_value = "2G")]
    mem: String,
    /// Directory for the spill file (default: the OS temp dir)
    #[arg(long)]
    temp: Option<PathBuf>,
    #[clap(flatten)]
    verbose: clap_verbosity_flag::Verbosity,
}

fn main() -> io::Result<()> {
    let cli = Cli::parse();
    env_logger::Builder::new()
        .filter_level(cli.verbose.log_level_filter())
        .init();
    let t0 = Instant::now();

    let mut dem = Dem::new(
        cli.res,
        parse_origin(&cli.origin)?,
        cli.tile_size,
        cli.max_edge,
        parse_mem(&cli.mem)?,
        cli.temp.clone(),
        cli.bbox.as_deref().map(parse_bbox).transpose()?,
    );
    //-- fail fast on a bad --bbox before any streaming
    if dem.user_bbox.is_some() {
        dem.ensure_grid()?;
    }

    let stdin = io::stdin();
    let mut r = BufReader::with_capacity(1 << 20, stdin.lock());
    if sst::stream::peek_is_binary(&mut r, sst::stream::SMA_MAGIC)? {
        info!("input stream: binary .sma");
        sst::stream::skip_magic(&mut r, sst::stream::SMA_MAGIC)?;
        while let Some(rec) = sst::stream::sma_read(&mut r)? {
            dem.apply(rec)?;
        }
    } else {
        info!("input stream: ASCII .sma");
        for line in r.lines() {
            let l = line?;
            if l.is_empty() || l.starts_with('#') {
                continue;
            }
            match parse_sma_line(&l) {
                Some(rec) => dem.apply(rec)?,
                None => {
                    error!("Wrongly formatted stream. Abort.");
                    std::process::exit(1);
                }
            }
        }
    }

    dem.finish(&cli, t0)
}

/// The DEM being built: the vertex table of the streaming TIN, the grid and
/// its tile store, and the running stats.
struct Dem {
    //-- configuration
    res: f64,
    origin: (f64, f64),
    tile: u32,
    max_edge: f64,
    mem_budget: u64,
    temp_dir: Option<PathBuf>,
    user_bbox: Option<[f64; 4]>,
    //-- state
    vs: HashMap<u32, [f64; 3]>,
    next_id: u32,
    stream_bbox: Option<[f64; 4]>,
    stream_crs: (u32, u32),
    grid: Option<Grid>,
    store: Option<TileStore>,
    stats: RasterStats,
    vertices_in: u64,
    max_active: usize,
}

impl Dem {
    #[allow(clippy::too_many_arguments)]
    fn new(
        res: f64,
        origin: (f64, f64),
        tile: u32,
        max_edge: f64,
        mem_budget: u64,
        temp_dir: Option<PathBuf>,
        user_bbox: Option<[f64; 4]>,
    ) -> Dem {
        Dem {
            res,
            origin,
            tile,
            max_edge,
            mem_budget,
            temp_dir,
            user_bbox,
            vs: HashMap::new(),
            next_id: 1,
            stream_bbox: None,
            stream_crs: (0, 0),
            grid: None,
            store: None,
            stats: RasterStats::default(),
            vertices_in: 0,
            max_active: 0,
        }
    }

    fn apply(&mut self, rec: SmaRecord) -> io::Result<()> {
        match rec {
            SmaRecord::B(b) => {
                if self.stream_bbox.is_none() {
                    self.stream_bbox = Some(b);
                }
            }
            SmaRecord::P(h, v) => self.stream_crs = (h, v),
            SmaRecord::V(x, y, z) => {
                //-- ids are 1-based (ttdt numbers them from 1)
                self.vs.insert(self.next_id, [x, y, z]);
                self.next_id += 1;
                self.vertices_in += 1;
                self.max_active = self.max_active.max(self.vs.len());
            }
            SmaRecord::X(id) => {
                //-- the vertex will not be referenced by any future triangle
                self.vs.remove(&id);
            }
            SmaRecord::F(a, b, c) => {
                self.ensure_grid()?;
                let v0 = *self
                    .vs
                    .get(&a)
                    .ok_or_else(|| unknown_vertex(a))?;
                let v1 = *self
                    .vs
                    .get(&b)
                    .ok_or_else(|| unknown_vertex(b))?;
                let v2 = *self
                    .vs
                    .get(&c)
                    .ok_or_else(|| unknown_vertex(c))?;
                let grid = self.grid.as_ref().unwrap();
                let store = self.store.as_mut().unwrap();
                raster::rasterize(grid, store, [v0, v1, v2], self.max_edge, &mut self.stats)?;
            }
        }
        Ok(())
    }

    /// Create the grid and its tile store at the first triangle (the bbox and
    /// the CRS are only known once the stream header has been read).
    fn ensure_grid(&mut self) -> io::Result<()> {
        if self.grid.is_some() {
            return Ok(());
        }
        let bbox = self
            .user_bbox
            .or(self.stream_bbox)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "no bbox: the .sma stream has no `b` record and no --bbox was given",
                )
            })?;
        let grid = Grid::new(bbox, self.res, self.origin, self.tile)?;
        info!(
            "grid: {} x {} cells (res {}; bbox {:?})",
            grid.ncols,
            grid.nrows,
            grid.res,
            [grid.minx, grid.miny, grid.maxx, grid.maxy]
        );
        self.store = Some(TileStore::new(
            &grid,
            self.mem_budget,
            self.temp_dir.as_deref(),
            NODATA,
        )?);
        self.grid = Some(grid);
        Ok(())
    }

    fn finish(mut self, cli: &Cli, t0: Instant) -> io::Result<()> {
        self.ensure_grid()?;
        let grid = self.grid.unwrap();
        let mut store = self.store.unwrap();

        //-- CRS: the --crs flag overrides the stream's `p` record
        let (mut h, mut v) = self.stream_crs;
        let mut crs_source = "stream";
        if let Some(s) = &cli.crs {
            h = parse_epsg(s)? as u32;
            v = 0;
            crs_source = "--crs";
        }
        if let Some(s) = &cli.vertical_crs {
            v = parse_epsg(s)? as u32;
            crs_source = "flag(s)";
        }
        let epsg = epsg16(h);
        let vertical_epsg = if epsg.is_some() { epsg16(v) } else { None };
        if epsg.is_none() {
            eprintln!(
                "warning: no CRS found in the stream and no --crs given; \
                 writing the GeoTIFF without georeferencing"
            );
        } else if (4000..=4999).contains(&epsg.unwrap()) {
            eprintln!(
                "warning: EPSG:{} looks like a geographic (degree) CRS; \
                 --res and --max-edge are then in degrees, and linear \
                 interpolation in lon/lat distorts distances",
                epsg.unwrap()
            );
        } else {
            info!("CRS: EPSG:{} (from {})", epsg.unwrap(), crs_source);
        }

        let (compression, predictor) = parse_compression(&cli.compress)?;
        let overview_levels = if cli.overviews == "none" {
            Vec::new()
        } else {
            writer::auto_overview_levels(&grid)
        };
        let spec = writer::CogSpec {
            path: cli.output.clone(),
            epsg,
            vertical_epsg,
            compression,
            predictor,
            overview_levels,
        };
        writer::write_cog(&spec, &grid, &mut store)?;

        //-- stats
        let sstats = store.stats();
        let coverage = self.stats.cells_written as f64 / grid.full_cells() as f64 * 100.0;
        let out_bytes = std::fs::metadata(&cli.output).map(|m| m.len()).unwrap_or(0);
        info!("vertices read: {} (max {} active)", self.vertices_in, self.max_active);
        info!(
            "triangles: {} read, {} rasterised, {} skipped ({} degenerate, {} long-edge, {} outside)",
            self.stats.triangles,
            self.stats.rasterised,
            self.stats.skipped_degenerate + self.stats.skipped_long_edge + self.stats.skipped_outside,
            self.stats.skipped_degenerate,
            self.stats.skipped_long_edge,
            self.stats.skipped_outside
        );
        info!("longest 2D edge seen: {:.3}", self.stats.max_edge);
        info!(
            "cells written: {} of {} ({:.1}%)",
            self.stats.cells_written,
            grid.full_cells(),
            coverage
        );
        info!(
            "tile store: {} ({} tiles touched, {} evictions, {} reloads, spill {:.1} MB)",
            sstats.mode,
            sstats.tiles_touched,
            sstats.evictions,
            sstats.reloads,
            sstats.spill_bytes as f64 / (1024.0 * 1024.0)
        );
        info!("output: {} ({:.1} MB)", cli.output.display(), out_bytes as f64 / (1024.0 * 1024.0));
        info!("done in {:.1?}", t0.elapsed());
        if self.stats.cells_written == 0 {
            eprintln!("warning: not a single cell was written; check --res, --bbox and --max-edge");
        }
        Ok(())
    }
}

fn unknown_vertex(id: u32) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("triangle references unknown/finalised vertex id {id}"),
    )
}

/// "28992" or "EPSG:28992" -> 28992 (0 is not a valid EPSG code).
fn parse_epsg(s: &str) -> io::Result<u16> {
    let t = s.trim().to_lowercase();
    let t = t.strip_prefix("epsg:").unwrap_or(&t);
    let code: u16 = t
        .parse()
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, format!("bad EPSG code: {s:?}")))?;
    if code == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "EPSG code 0 means 'unknown'; give a real code or leave --crs out",
        ));
    }
    Ok(code)
}

/// Stream EPSG codes can be any u32; GeoTIFF wants u16.
fn epsg16(code: u32) -> Option<u16> {
    u16::try_from(code).ok().filter(|&c| c != 0)
}

fn parse_bbox(v: &[f64]) -> io::Result<[f64; 4]> {
    if v.len() != 4 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("--bbox takes 4 values (minx miny maxx maxy), got {}", v.len()),
        ));
    }
    Ok([v[0], v[1], v[2], v[3]])
}

fn parse_origin(v: &[f64]) -> io::Result<(f64, f64)> {
    if v.len() != 2 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("--origin takes 2 values (x y), got {}", v.len()),
        ));
    }
    Ok((v[0], v[1]))
}

/// "2G" -> 2 GiB; also accepts M, K, or a plain byte count.
fn parse_mem(s: &str) -> io::Result<u64> {
    let s = s.trim();
    let (num, mult) = match s.chars().last() {
        Some(c) if c.is_ascii_alphabetic() => {
            let m = match c.to_ascii_lowercase() {
                'k' => 1024u64,
                'm' => 1024 * 1024,
                'g' => 1024 * 1024 * 1024,
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("bad memory size: {s:?} (use eg 2G, 512M, 65536K)"),
                    ))
                }
            };
            (&s[..s.len() - 1], m)
        }
        _ => (s, 1),
    };
    let n: u64 = num
        .trim()
        .parse()
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, format!("bad memory size: {s:?}")))?;
    Ok(n.checked_mul(mult).ok_or_else(|| {
        io::Error::new(io::ErrorKind::InvalidInput, format!("memory size {s:?} overflows"))
    })?)
}

fn parse_compression(s: &str) -> io::Result<(Compression, Predictor)> {
    match s.to_lowercase().as_str() {
        "deflate" => Ok((Compression::Deflate, Predictor::FloatingPoint)),
        "lzw" => Ok((Compression::Lzw, Predictor::FloatingPoint)),
        "none" => Ok((Compression::None, Predictor::None)),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("bad --compress: {s:?} (deflate | lzw | none)"),
        )),
    }
}

/// One `.sma` record from an ASCII line (the binary path decodes the same
/// enum directly).
fn parse_sma_line(l: &str) -> Option<SmaRecord> {
    let mut it = l.split_whitespace();
    match it.next()? {
        "b" => Some(SmaRecord::B([
            it.next()?.parse().ok()?,
            it.next()?.parse().ok()?,
            it.next()?.parse().ok()?,
            it.next()?.parse().ok()?,
        ])),
        "v" => Some(SmaRecord::V(
            it.next()?.parse().ok()?,
            it.next()?.parse().ok()?,
            it.next()?.parse().ok()?,
        )),
        "f" => Some(SmaRecord::F(
            it.next()?.parse().ok()?,
            it.next()?.parse().ok()?,
            it.next()?.parse().ok()?,
        )),
        "x" => Some(SmaRecord::X(it.next()?.parse().ok()?)),
        "p" => Some(SmaRecord::P(
            it.next()?.parse().ok()?,
            it.next()?.parse().ok()?,
        )),
        _ => None,
    }
}
