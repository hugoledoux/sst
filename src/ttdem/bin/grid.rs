//-- output grid: georeferencing + cell/tile indexing
//--
//-- Rows go north-to-south (row 0 = north, GeoTIFF convention); pixel
//-- *centres* are at (minx + (c+0.5)*res, maxy - (r+0.5)*res), and the
//-- raster bbox is the outer edge of the pixel areas (PixelIsArea).

use std::io;

//-- slack (in cell units) when rounding span boundaries, so that samples
//-- that land exactly on a triangle edge are never lost to fp noise.
//-- Duplicate writes on shared edges are harmless (both planes agree there).
pub const SPAN_EPS: f64 = 1e-6;

#[derive(Clone, Copy, Debug)]
pub struct Grid {
    pub minx: f64,
    pub miny: f64,
    pub maxx: f64,
    pub maxy: f64,
    pub res: f64,
    pub ncols: u32,
    pub nrows: u32,
    pub tile: u32,
    pub tiles_x: u32,
    pub tiles_y: u32,
}

impl Grid {
    /// Snap a bbox outward onto the grid aligned at `origin` (the default
    /// origin is (0, 0) in CRS units, which keeps neighbouring runs aligned).
    pub fn new(bbox: [f64; 4], res: f64, origin: (f64, f64), tile: u32) -> io::Result<Grid> {
        if !res.is_finite() || res <= 0.0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("--res must be a positive number, got {res}"),
            ));
        }
        if !bbox.iter().all(|v| v.is_finite()) || bbox[0] >= bbox[2] || bbox[1] >= bbox[3] {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("invalid bbox [{:?}] (need minx < maxx and miny < maxy)", bbox),
            ));
        }
        let snap_lo = |v: f64, o: f64| o + ((v - o) / res).floor() * res;
        let snap_hi = |v: f64, o: f64| o + ((v - o) / res).ceil() * res;
        let minx = snap_lo(bbox[0], origin.0);
        let mut maxx = snap_hi(bbox[2], origin.0);
        let miny = snap_lo(bbox[1], origin.1);
        let mut maxy = snap_hi(bbox[3], origin.1);
        if !minx.is_finite() || !maxy.is_finite() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "bbox/resolution combination overflows",
            ));
        }
        let ncols = (((maxx - minx) / res).round() as u64).max(1);
        let nrows = (((maxy - miny) / res).round() as u64).max(1);
        if ncols > u32::MAX as u64 || nrows > u32::MAX as u64 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("grid {ncols} x {nrows} exceeds the GeoTIFF u32 dimensions"),
            ));
        }
        //-- recompute the corners so they are exact multiples of res
        maxx = minx + (ncols as f64) * res;
        maxy = miny + (nrows as f64) * res;
        let tiles_x = (ncols as u32).div_ceil(tile);
        let tiles_y = (nrows as u32).div_ceil(tile);
        Ok(Grid {
            minx,
            miny,
            maxx,
            maxy,
            res,
            ncols: ncols as u32,
            nrows: nrows as u32,
            tile,
            tiles_x,
            tiles_y,
        })
    }

    pub fn center_x(&self, c: u32) -> f64 {
        self.minx + (c as f64 + 0.5) * self.res
    }

    pub fn center_y(&self, r: u32) -> f64 {
        self.maxy - (r as f64 + 0.5) * self.res
    }

    /// Inclusive column range whose cell centres lie in [xmin, xmax], clamped
    /// to the grid; `None` if the range is empty.
    pub fn col_span(&self, xmin: f64, xmax: f64) -> Option<(u32, u32)> {
        let a = ((xmin - self.minx) / self.res - 0.5 - SPAN_EPS).ceil() as i64;
        let b = ((xmax - self.minx) / self.res - 0.5 + SPAN_EPS).floor() as i64;
        span_from(a, b, self.ncols as i64)
    }

    /// Inclusive row range (north-to-south) whose cell centres lie in
    /// [ymin, ymax], clamped to the grid; `None` if the range is empty.
    pub fn row_span(&self, ymin: f64, ymax: f64) -> Option<(u32, u32)> {
        let a = ((self.maxy - ymax) / self.res - 0.5 - SPAN_EPS).ceil() as i64;
        let b = ((self.maxy - ymin) / self.res - 0.5 + SPAN_EPS).floor() as i64;
        span_from(a, b, self.nrows as i64)
    }

    /// Size in bytes of the full-resolution raster; `None` on overflow
    /// (which simply forces the spilled store).
    pub fn full_bytes(&self) -> Option<u64> {
        self.full_cells().checked_mul(4)
    }

    pub fn full_cells(&self) -> u64 {
        self.ncols as u64 * self.nrows as u64
    }

    pub fn tile_cells(&self) -> usize {
        self.tile as usize * self.tile as usize
    }
}

fn span_from(a: i64, b: i64, limit: i64) -> Option<(u32, u32)> {
    //-- a span entirely outside the raster is not a span
    if b < 0 || a >= limit {
        return None;
    }
    let a = a.clamp(0, limit - 1);
    let b = b.clamp(0, limit - 1);
    if a > b {
        None
    } else {
        Some((a as u32, b as u32))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn snapping_is_outward_and_origin_aligned() {
        //-- 1km-ish at 50cm from an off-origin bbox
        let g = Grid::new([123.4, 456.7, 1123.9, 1456.2], 0.5, (0.0, 0.0), 512).unwrap();
        assert_eq!(g.minx, 123.0);
        assert_eq!(g.maxx, 1124.0);
        assert_eq!(g.ncols, 2002);
        assert_eq!((g.maxx - g.minx) / g.res, 2002.0);
        assert_eq!(g.tiles_x, 4);
        assert_eq!(g.nrows, 2000);
    }

    #[test]
    fn spans_hit_cells_whose_centre_is_inside() {
        let g = Grid::new([0.0, 0.0, 10.0, 10.0], 1.0, (0.0, 0.0), 512).unwrap();
        //-- centres at 0.5..9.5
        assert_eq!(g.col_span(0.0, 10.0), Some((0, 9)));
        assert_eq!(g.col_span(2.4, 3.6), Some((2, 3)));
        assert_eq!(g.col_span(2.6, 3.4), None); //-- no centre in between
        assert_eq!(g.col_span(2.6, 2.6), None);
        assert_eq!(g.col_span(8.4, 12.0), Some((8, 9))); //-- clamped at the edge
        assert_eq!(g.row_span(0.0, 10.0), Some((0, 9)));
    }

    #[test]
    fn tiling_covers_edges() {
        let g = Grid::new([0.0, 0.0, 1000.0, 1000.0], 1.0, (0.0, 0.0), 512).unwrap();
        assert_eq!((g.tiles_x, g.tiles_y), (2, 2));
        assert_eq!(g.tile_cells(), 512 * 512);
    }
}
