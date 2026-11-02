//-- linear interpolation of one TIN triangle onto the raster grid
//--
//-- For each triangle: fit the plane z(x, y) through its three vertices,
//-- then rasterise it scanline by scanline -- for every raster row that
//-- crosses the triangle, the covered cell centres form one contiguous
//-- span, so no per-cell point-in-triangle test is needed. Elevation is
//-- evaluated once at the start of the span and advanced by a*res per
//-- column.

use std::io;

use crate::grid::Grid;
use crate::store::TileStore;

#[derive(Default)]
pub struct RasterStats {
    pub triangles: u64,
    pub rasterised: u64,
    pub skipped_degenerate: u64,
    pub skipped_long_edge: u64,
    pub skipped_outside: u64,
    pub cells_written: u64,
    pub max_edge: f64,
}

/// Rasterise one triangle (as three xyz vertices, any winding).
pub fn rasterize(
    grid: &Grid,
    store: &mut TileStore,
    tri: [[f64; 3]; 3],
    max_edge: f64,
    stats: &mut RasterStats,
) -> io::Result<()> {
    stats.triangles += 1;

    //-- reject triangles that bridge data gaps (Isenburg §3.1); 0 disables.
    //-- the longest edge seen is tracked regardless, so a threshold can be
    //-- picked from the data
    for i in 0..3 {
        let p = tri[i];
        let q = tri[(i + 1) % 3];
        let e = ((q[0] - p[0]).powi(2) + (q[1] - p[1]).powi(2)).sqrt();
        if e > stats.max_edge {
            stats.max_edge = e;
        }
        if max_edge > 0.0 && e > max_edge {
            stats.skipped_long_edge += 1;
            return Ok(());
        }
    }

    let [v0, v1, v2] = tri;
    let dx1 = v1[0] - v0[0];
    let dy1 = v1[1] - v0[1];
    let dz1 = v1[2] - v0[2];
    let dx2 = v2[0] - v0[0];
    let dy2 = v2[1] - v0[1];
    let dz2 = v2[2] - v0[2];

    //-- twice the signed area: skip exactly-collinear triangles and any
    //-- with a non-finite z. Thin slivers are kept: inside them the plane
    //-- still interpolates exactly between the three vertices.
    let det = dx1 * dy2 - dx2 * dy1;
    if det == 0.0 || tri.iter().any(|v| !v[2].is_finite()) {
        stats.skipped_degenerate += 1;
        return Ok(());
    }
    let a = (dz1 * dy2 - dz2 * dy1) / det; //-- dz/dx
    let b = (dz2 * dx1 - dz1 * dx2) / det; //-- dz/dy
    if !a.is_finite() || !b.is_finite() {
        stats.skipped_degenerate += 1;
        return Ok(());
    }

    let (r0, r1) = match grid.row_span(
        v0[1].min(v1[1]).min(v2[1]),
        v0[1].max(v1[1]).max(v2[1]),
    ) {
        Some(s) => s,
        None => {
            stats.skipped_outside += 1;
            return Ok(());
        }
    };

    for r in r0..=r1 {
        let y = grid.center_y(r);
        //-- intersect the horizontal scanline with the 3 edges; the crossings
        //-- of a convex polygon give the covered x-span
        let mut xmin = f64::INFINITY;
        let mut xmax = f64::NEG_INFINITY;
        for i in 0..3 {
            let p = tri[i];
            let q = tri[(i + 1) % 3];
            if (p[1] - y) * (q[1] - y) > 0.0 {
                continue; //-- edge entirely above or below the scanline
            }
            if p[1] == q[1] {
                //-- horizontal edge lying on the scanline: both endpoints
                if p[1] == y {
                    xmin = xmin.min(p[0]).min(q[0]);
                    xmax = xmax.max(p[0]).max(q[0]);
                }
                continue;
            }
            let t = (y - p[1]) / (q[1] - p[1]);
            if !t.is_finite() {
                continue;
            }
            let x = p[0] + t * (q[0] - p[0]);
            if x < xmin {
                xmin = x;
            }
            if x > xmax {
                xmax = x;
            }
        }
        if xmin > xmax {
            continue;
        }
        let (c0, c1) = match grid.col_span(xmin, xmax) {
            Some(s) => s,
            None => continue,
        };

        //-- evaluate z once, then step a*res per column
        let x0 = grid.center_x(c0);
        let z_start = v0[2] + a * (x0 - v0[0]) + b * (y - v0[1]);
        let dz = a * grid.res;
        store.span_write(grid, r, c0, c1, z_start, dz)?;
        stats.cells_written += (c1 - c0 + 1) as u64;
    }
    stats.rasterised += 1;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(tris: Vec<[[f64; 3]; 3]>) -> (Vec<f32>, Grid, RasterStats) {
        let grid = Grid::new([0.0, 0.0, 10.0, 10.0], 1.0, (0.0, 0.0), 512).unwrap();
        let mut store = TileStore::new_full(&grid, -9999.0);
        let mut stats = RasterStats::default();
        for t in tris {
            rasterize(&grid, &mut store, t, 0.0, &mut stats).unwrap();
        }
        (store.full_grid(), grid, stats)
    }

    #[test]
    fn two_triangles_tile_a_10x10_plane() {
        //-- the unit square split along the diagonal: z = max(x, y)
        let (data, grid, stats) = run(vec![
            [[0.0, 0.0, 0.0], [10.0, 0.0, 10.0], [10.0, 10.0, 10.0]],
            [[0.0, 0.0, 0.0], [10.0, 10.0, 10.0], [0.0, 10.0, 10.0]],
        ]);
        assert_eq!(stats.rasterised, 2);
        //-- every cell covered; on-diagonal cells are written by both
        //-- triangles (same z there), so the count is >= 100
        assert!(stats.cells_written >= 100);
        for r in 0..grid.nrows as usize {
            for c in 0..grid.ncols as usize {
                let x = grid.center_x(c as u32);
                let y = grid.center_y(r as u32);
                assert_eq!(
                    data[r * grid.ncols as usize + c], x.max(y) as f32,
                    "cell ({r},{c})"
                );
            }
        }
    }

    #[test]
    fn every_cell_is_covered_exactly_across_shared_edges() {
        //-- two triangles tiling the square [0,4]²; no hole may be left
        //-- inside. With row 0 = north, the diamond sits in the
        //-- bottom-left of the raster
        let (data, grid, stats) = run(vec![
            [[0.0, 0.0, 1.0], [4.0, 0.0, 1.0], [0.0, 4.0, 1.0]],
            [[4.0, 4.0, 1.0], [0.0, 4.0, 1.0], [4.0, 0.0, 1.0]],
        ]);
        assert_eq!(stats.rasterised, 2);
        let ncols = grid.ncols as usize;
        //-- the diamond covers world y in [0,4] = raster rows 6..9
        for r in 6..10usize {
            for c in 0..4usize {
                assert_eq!(data[r * ncols + c], 1.0, "hole at cell ({r},{c})");
            }
        }
        //-- and the far corner of the raster is untouched nodata
        assert_eq!(data[ncols + 0], -9999.0); //-- row 1, col 0: north-west
        assert_eq!(data[0], -9999.0);
    }

    #[test]
    fn outside_and_degenerate_are_counted() {
        let (_, _g, stats) = run(vec![
            [[100.0, 100.0, 0.0], [101.0, 100.0, 0.0], [100.0, 101.0, 0.0]],
            [[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 2.0, 0.0]],
        ]);
        assert_eq!(stats.skipped_outside, 1);
        assert_eq!(stats.skipped_degenerate, 1);
    }

    #[test]
    fn long_edges_are_rejected_when_enabled() {
        let grid = Grid::new([0.0, 0.0, 10.0, 10.0], 1.0, (0.0, 0.0), 512).unwrap();
        let mut store = TileStore::new_full(&grid, -9999.0);
        let mut stats = RasterStats::default();
        rasterize(
            &grid,
            &mut store,
            [[0.0, 0.0, 0.0], [10.0, 0.0, 5.0], [10.0, 10.0, 5.0]],
            5.0,
            &mut stats,
        )
        .unwrap();
        assert_eq!(stats.skipped_long_edge, 1);
        assert_eq!(stats.rasterised, 0);
    }
}
