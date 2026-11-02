//-- COG output for ttdem, on top of the geotiff-writer crate
//--
//-- The tile store hands over complete tile-row-major tiles; this module
//-- wraps them in the COG container (BigTIFF, deflate + floating-point
//-- predictor, georeferencing, internal overviews). The output is written to
//-- `<path>.tmp` and renamed on success, so a failed run never leaves a
//-- half-written COG behind.

use std::io::{self, BufWriter, Write as _};
use std::path::PathBuf;

use geotiff_writer::cog::CogBuilder;
use geotiff_writer::{Compression, GeoTiffBuilder, Predictor, RasterType};

pub struct CogSpec {
    pub path: PathBuf,
    /// horizontal + vertical EPSG; `None` skips georeferencing
    pub epsg: Option<u16>,
    pub vertical_epsg: Option<u16>,
    pub compression: Compression,
    pub predictor: Predictor,
    /// overview subsampling factors (powers of two); empty = no overviews
    pub overview_levels: Vec<u32>,
}

/// Generate the overview levels for a grid: halve the largest dimension
/// until it is no bigger than a tile.
pub fn auto_overview_levels(grid: &crate::grid::Grid) -> Vec<u32> {
    let max_dim = grid.ncols.max(grid.nrows) as u64;
    let tile = grid.tile as u64;
    let mut levels = Vec::new();
    let mut factor = 2u64;
    while max_dim / factor > tile && levels.len() < 12 {
        levels.push(factor as u32);
        factor *= 2;
    }
    levels
}

pub fn write_cog(
    spec: &CogSpec,
    grid: &crate::grid::Grid,
    store: &mut crate::store::TileStore,
) -> io::Result<()> {
    let tile = grid.tile as usize;
    let tiff = GeoTiffBuilder::new(grid.ncols, grid.nrows)
        .tile_size(grid.tile, grid.tile)
        .compression(spec.compression)
        .predictor(spec.predictor)
        .raster_type(RasterType::PixelIsArea);
    let tiff = match spec.epsg {
        Some(e) => tiff.epsg(e),
        None => tiff,
    };
    let tiff = match spec.vertical_epsg {
        Some(v) => tiff.vertical_epsg(v),
        None => tiff,
    };
    let cog = if spec.overview_levels.is_empty() {
        CogBuilder::new(tiff).no_overviews()
    } else {
        CogBuilder::new(tiff).overview_levels(spec.overview_levels.clone())
    };

    //-- write into a temp file next to the output, then rename
    let dir = spec
        .path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| PathBuf::from("."));
    let tmp = tempfile::Builder::new()
        .prefix("ttdem-")
        .suffix(".tmp")
        .tempfile_in(&dir)?;

    {
        //-- the tile buffers are strided sub-rectangles (the full-grid store
        //-- hands out slices of the whole raster, edge tiles are partial), so
        //-- compact each tile into a contiguous buffer before writing; the
        //-- copy is negligible next to the compression that follows
        let sink = BufWriter::new(tmp.reopen()?);
        let mut tw = cog.tile_writer(sink).map_err(map_err)?;
        let (ncols, nrows) = (grid.ncols as usize, grid.nrows as usize);
        let mut scratch: Vec<f32> = Vec::new();
        store.for_each_tile(grid, &mut |tr, tc, data, stride| {
            let h = (nrows - tr as usize * tile).min(tile);
            let w = (ncols - tc as usize * tile).min(tile);
            scratch.clear();
            for r in 0..h {
                scratch.extend_from_slice(&data[r * stride..r * stride + w]);
            }
            let view = ndarray::ArrayView2::from_shape((h, w), &scratch)
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;
            tw.write_tile(tc as usize * tile, tr as usize * tile, &view)
                .map_err(map_err)?;
            Ok(())
        })?;
        let mut sink = tw.finish().map_err(map_err)?;
        sink.flush()?; //-- then drop; the file moves into place below
    }

    tmp.persist(&spec.path).map_err(|e| {
        io::Error::new(
            io::ErrorKind::Other,
            format!("could not move the temp COG into place: {}", e.error),
        )
    })?;
    Ok(())
}

fn map_err(e: geotiff_writer::Error) -> io::Error {
    io::Error::new(io::ErrorKind::Other, format!("COG writer: {}", e))
}
