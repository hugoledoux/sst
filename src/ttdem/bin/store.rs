//-- tile store: where the rasterised cells live while the TIN streams in
//--
//-- Two backends behind one interface:
//--   * Full    -- the whole raster as one Vec<f32>; used when it fits the
//--               --mem budget. Simple, fast, zero-copy tile handover.
//--   * Spilled -- an LRU cache of tile buffers; evicted tiles are
//--               deflate-compressed and appended to one spill file (with an
//--               index), and read back when a later triangle touches them
//--               again. This is the 2026 analogue of Isenburg's row-band
//--               temp files: the tile is the bucket, no global sort needed.
//--
//-- Both backends hand tiles to the COG writer at the end (in the same order,
//-- tile row-major); the writer does the compression/IFD work. The two
//-- backends must produce byte-identical COGs for the same input (verified by
//-- feeding the same .sma twice with different --mem).

use std::collections::BTreeMap;
use std::io::{self, Read, Seek, SeekFrom, Write};

use flate2::read::DeflateDecoder;
use flate2::write::DeflateEncoder;
use flate2::Compression as DeflateLevel;
use hashbrown::HashMap;

use crate::grid::Grid;

pub struct StoreStats {
    pub mode: &'static str,
    pub evictions: u64,
    pub reloads: u64,
    pub spill_bytes: u64,
    pub tiles_touched: u64,
}

pub enum TileStore {
    Full(FullStore),
    Spilled(Box<SpillStore>),
}

impl TileStore {
    #[cfg(test)]
    /// Force the full-RAM backend, regardless of --mem (used by the tests).
    pub fn new_full(grid: &Grid, nodata: f32) -> TileStore {
        TileStore::Full(FullStore::new(grid, nodata))
    }

    /// Full store when the whole raster fits the memory budget, spilled
    /// tile cache otherwise.
    pub fn new(
        grid: &Grid,
        mem_budget: u64,
        temp_dir: Option<&std::path::Path>,
        nodata: f32,
    ) -> io::Result<TileStore> {
        match grid.full_bytes() {
            Some(bytes) if bytes <= mem_budget => {
                Ok(TileStore::Full(FullStore::new(grid, nodata)))
            }
            _ => Ok(TileStore::Spilled(Box::new(SpillStore::new(
                grid, mem_budget, temp_dir, nodata,
            )?))),
        }
    }

    /// Write one horizontal run of cells [c0, c1] on `row`; `z_start` is the
    /// interpolated elevation at c0 and `dz` the step per column.
    pub fn span_write(
        &mut self,
        grid: &Grid,
        row: u32,
        c0: u32,
        c1: u32,
        z_start: f64,
        dz: f64,
    ) -> io::Result<()> {
        match self {
            TileStore::Full(s) => Ok(s.span_write(grid, row, c0, c1, z_start, dz)),
            TileStore::Spilled(s) => s.span_write(grid, row, c0, c1, z_start, dz),
        }
    }

    /// Hand every tile over to the callback exactly once, in tile-row-major
    /// order. `data[0]` is the tile's first sample and `stride` the number of
    /// samples per row, so the caller can build a (h, w) view with strides
    /// (stride, 1); h and w are derived by the caller from the grid. The full
    /// store hands out *all* tiles (untouched ones are nodata), the spilled
    /// store only the touched ones (the writer fills the rest with nodata).
    pub fn for_each_tile(
        &mut self,
        grid: &Grid,
        cb: &mut dyn FnMut(u32, u32, &[f32], usize) -> io::Result<()>,
    ) -> io::Result<()> {
        match self {
            TileStore::Full(s) => s.for_each_tile(grid, cb),
            TileStore::Spilled(s) => s.for_each_tile(grid, cb),
        }
    }

    pub fn stats(&self) -> StoreStats {
        match self {
            TileStore::Full(_) => StoreStats {
                mode: "full",
                evictions: 0,
                reloads: 0,
                spill_bytes: 0,
                tiles_touched: 0,
            },
            TileStore::Spilled(s) => StoreStats {
                mode: "spilled",
                evictions: s.evictions,
                reloads: s.reloads,
                spill_bytes: s.spill_end,
                tiles_touched: s.index.len() as u64,
            },
        }
    }

    #[cfg(test)]
    pub fn full_grid(&mut self) -> Vec<f32> {
        match self {
            TileStore::Full(s) => s.data.clone(),
            _ => unreachable!(),
        }
    }
}

//-- full raster ----------------------------------------------------------------

pub struct FullStore {
    data: Vec<f32>,
}

impl FullStore {
    fn new(grid: &Grid, nodata: f32) -> FullStore {
        FullStore {
            data: vec![nodata; grid.full_cells() as usize],
        }
    }

    fn span_write(&mut self, grid: &Grid, row: u32, c0: u32, c1: u32, z_start: f64, dz: f64) {
        let base = row as usize * grid.ncols as usize;
        let mut z = z_start;
        for c in c0..=c1 {
            self.data[base + c as usize] = z as f32;
            z += dz;
        }
    }

    fn for_each_tile(
        &mut self,
        grid: &Grid,
        cb: &mut dyn FnMut(u32, u32, &[f32], usize) -> io::Result<()>,
    ) -> io::Result<()> {
        let tile = grid.tile as usize;
        let ncols = grid.ncols as usize;
        for tr in 0..grid.tiles_y {
            for tc in 0..grid.tiles_x {
                let base = tr as usize * tile * ncols + tc as usize * tile;
                cb(tr, tc, &self.data[base..], ncols)?;
            }
        }
        Ok(())
    }
}

//-- LRU + spill ------------------------------------------------------------------

pub struct SpillStore {
    tile: u32,
    tile_cells: usize,
    tile_bytes: usize,
    budget: usize,
    nodata: f32,
    resident: HashMap<u64, (Box<[f32]>, u64)>,
    order: BTreeMap<u64, u64>,
    clock: u64,
    resident_bytes: usize,
    spill: tempfile::NamedTempFile,
    index: HashMap<u64, (u64, u32)>,
    spill_end: u64,
    evictions: u64,
    reloads: u64,
}

impl SpillStore {
    fn new(
        grid: &Grid,
        mem_budget: u64,
        temp_dir: Option<&std::path::Path>,
        nodata: f32,
    ) -> io::Result<SpillStore> {
        let tile_cells = grid.tile_cells();
        let tile_bytes = tile_cells * 4;
        let mut builder = tempfile::Builder::new();
        builder.prefix("ttdem-spill-");
        let spill = match temp_dir {
            Some(d) => builder.tempfile_in(d)?,
            None => builder.tempfile()?,
        };
        Ok(SpillStore {
            tile: grid.tile,
            tile_cells,
            tile_bytes,
            budget: (mem_budget as usize).max(tile_bytes),
            nodata,
            resident: HashMap::new(),
            order: BTreeMap::new(),
            clock: 0,
            resident_bytes: 0,
            spill,
            index: HashMap::new(),
            spill_end: 0,
            evictions: 0,
            reloads: 0,
        })
    }

    fn touch(&mut self, key: u64) {
        self.clock += 1;
        if let Some(entry) = self.resident.get_mut(&key) {
            let old = entry.1;
            entry.1 = self.clock;
            self.order.remove(&old);
            self.order.insert(self.clock, key);
        }
    }

    /// The buffer of one tile, loading it back from the spill file (or
    /// allocating it fresh) when it is not resident.
    fn tile_mut(&mut self, key: u64) -> io::Result<&mut Box<[f32]>> {
        if !self.resident.contains_key(&key) {
            let buf: Vec<f32> = match self.index.get(&key) {
                Some(&(off, len)) => {
                    let mut raw = vec![0u8; len as usize];
                    let f = self.spill.as_file_mut();
                    f.seek(SeekFrom::Start(off))?;
                    f.read_exact(&mut raw)?;
                    decode_tile(&raw, self.tile_cells)?
                }
                None => vec![self.nodata; self.tile_cells],
            };
            self.resident.insert(key, (buf.into_boxed_slice(), 0));
            self.resident_bytes += self.tile_bytes;
            self.reloads += 1;
            self.touch(key);
            self.evict_if_needed()?;
        } else {
            self.touch(key);
        }
        Ok(&mut self.resident.get_mut(&key).unwrap().0)
    }

    fn evict_if_needed(&mut self) -> io::Result<()> {
        while self.resident_bytes > self.budget {
            //-- `order` maps clock -> tile key; the smallest clock is the LRU
            let Some((&_clock, &key)) = self.order.iter().next() else {
                break;
            };
            self.order.remove(&_clock);
            let (data, _) = self.resident.remove(&key).unwrap();
            self.resident_bytes -= self.tile_bytes;
            let blob = encode_tile(&data);
            match self.index.get(&key) {
                //-- rewrite in place when the new blob fits the old slot
                Some(&(off, len)) if blob.len() as u32 <= len => {
                    let f = self.spill.as_file_mut();
                    f.seek(SeekFrom::Start(off))?;
                    f.write_all(&blob)?;
                    self.index.insert(key, (off, blob.len() as u32));
                }
                _ => {
                    let off = self.spill_end;
                    let f = self.spill.as_file_mut();
                    f.seek(SeekFrom::Start(off))?;
                    f.write_all(&blob)?;
                    self.spill_end += blob.len() as u64;
                    self.index.insert(key, (off, blob.len() as u32));
                }
            }
            self.evictions += 1;
        }
        Ok(())
    }

    fn span_write(
        &mut self,
        grid: &Grid,
        row: u32,
        c0: u32,
        c1: u32,
        z_start: f64,
        dz: f64,
    ) -> io::Result<()> {
        let tile = grid.tile;
        let tr = row / tile;
        let lr = (row % tile) as usize;
        let mut z = z_start;
        let mut c = c0;
        while c <= c1 {
            let tc = c / tile;
            let seg_end = ((tc + 1) * tile - 1).min(c1);
            let key = ((tr as u64) << 32) | tc as u64;
            let buf = self.tile_mut(key)?;
            let t0 = (tc * tile) as usize;
            while c <= seg_end {
                buf[lr * tile as usize + (c as usize - t0)] = z as f32;
                z += dz;
                c += 1;
            }
        }
        Ok(())
    }

    fn for_each_tile(
        &mut self,
        _grid: &Grid,
        cb: &mut dyn FnMut(u32, u32, &[f32], usize) -> io::Result<()>,
    ) -> io::Result<()> {
        //-- resident first (taken out to avoid clones), then what is on disk;
        //-- sorted for a deterministic callback order
        let mut keys: Vec<u64> = self.resident.keys().copied().collect();
        for k in self.index.keys() {
            if !self.resident.contains_key(k) {
                keys.push(*k);
            }
        }
        keys.sort_unstable();
        let stride = self.tile as usize;
        for key in keys {
            let buf: Box<[f32]> = match self.resident.remove(&key) {
                Some((d, last)) => {
                    self.order.remove(&last);
                    self.resident_bytes -= self.tile_bytes;
                    d
                }
                None => {
                    let &(off, len) = self.index.get(&key).unwrap();
                    let mut raw = vec![0u8; len as usize];
                    let f = self.spill.as_file_mut();
                    f.seek(SeekFrom::Start(off))?;
                    f.read_exact(&mut raw)?;
                    decode_tile(&raw, self.tile_cells)?.into_boxed_slice()
                }
            };
            let tr = (key >> 32) as u32;
            let tc = (key & 0xffff_ffff) as u32;
            cb(tr, tc, &buf, stride)?;
        }
        Ok(())
    }
}

fn encode_tile(tile: &[f32]) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(tile.len() * 4);
    for v in tile {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    let mut enc = DeflateEncoder::new(Vec::new(), DeflateLevel::fast());
    let _ = enc.write_all(&bytes);
    enc.finish().unwrap()
}

fn decode_tile(raw: &[u8], tile_cells: usize) -> io::Result<Vec<f32>> {
    let mut bytes = vec![0u8; tile_cells * 4];
    DeflateDecoder::new(raw).read_exact(&mut bytes)?;
    Ok(bytes
        .chunks_exact(4)
        .map(|b| f32::from_le_bytes([b[0], b[1], b[2], b[3]]))
        .collect())
}
