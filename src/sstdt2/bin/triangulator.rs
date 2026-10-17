mod geom;

use std::fmt;
use std::fs::File;
use std::io::Write;
use std::io::{self};

use geojson::{Feature, FeatureCollection, Geometry, Value};
use serde_json::{to_value, Map};

use hashbrown::HashMap;
use std::collections::HashSet;

use rand::{thread_rng, Rng};

//-- number of triangles kept in the "most recently created" ring buffer;
//-- the walk's brute-force fallback scans these first, since the triangle
//-- containing a query point is almost always a recent one (Isenburg et al. §4.1)
const RECENT_TS_CAP: usize = 4096;

//----------------------
struct Qtcell {
    pts: Vec<usize>,
    ts: HashMap<usize, (usize, usize, usize)>,
    finalised: bool,
    walk_hints: Vec<usize>,
}

impl Qtcell {
    pub fn new() -> Qtcell {
        let p: Vec<usize> = Vec::new();
        let t: HashMap<usize, (usize, usize, usize)> = HashMap::new();
        let wh: Vec<usize> = Vec::new();
        Qtcell {
            pts: p,
            ts: t,
            finalised: false,
            walk_hints: wh,
        }
    }

    pub fn number_pts(&self) -> usize {
        self.pts.len()
    }

    pub fn number_ts(&self) -> usize {
        self.ts.len()
    }

    pub fn add_pt(&mut self, vi: usize) {
        self.pts.push(vi);
    }

    pub fn add_ts(&mut self, ti: usize, v0: usize, v1: usize, v2: usize) {
        self.ts.insert(ti, (v0, v1, v2));
    }

    pub fn set_final(&mut self, b: bool) {
        self.finalised = b
    }

    //-- keep a list of the most recently created triangles of this cell;
    //-- used to restart the walk when it falls into finalised space
    //-- (the mechanism prescribed by Isenburg et al. §4.1)
    pub fn add_walk_hint(&mut self, ti: usize) {
        if self.walk_hints.len() >= 16 {
            self.walk_hints.remove(0);
        }
        self.walk_hints.push(ti);
    }
}

//----------------------
// ┌────┬────┐
// │  1 │  3 │
// ├────┼────┤
// │  0 │  2 │
// └────┴────┘
//----------------------
struct Quadtree {
    cells: HashMap<Vec<u8>, Qtcell>,
    pub minx: f64,
    pub miny: f64,
    pub maxx: f64,
    pub maxy: f64,
    pub cellsize: usize,
    pub griddim: usize,
    depth: u32,
}

impl Quadtree {
    pub fn new() -> Quadtree {
        let cs: HashMap<Vec<u8>, Qtcell> = HashMap::new();
        Quadtree {
            cells: cs,
            cellsize: 0,
            griddim: 0,
            minx: std::f64::MAX,
            miny: std::f64::MAX,
            maxx: std::f64::MIN,
            maxy: std::f64::MIN,
            depth: 0,
        }
    }

    pub fn init(&mut self, griddim: usize) {
        self.griddim = griddim;
        self.depth = (griddim as f64).log(2.0).ceil() as u32;
        for i in 0..griddim {
            for j in 0..griddim {
                let qtc = self.get_qtc_from_gxgy(i, j);
                let nc = Qtcell::new();
                self.cells.insert(qtc, nc);
            }
        }
    }

    pub fn insert_one_vi(&mut self, x: f64, y: f64, vi: usize) {
        let q = self.get_qtc_from_xy(x, y);
        self.cells.get_mut(&q).unwrap().add_pt(vi);
    }

    pub fn get_cell_pts(&self, gx: usize, gy: usize) -> Vec<usize> {
        let q = self.get_qtc_from_gxgy(gx, gy);
        let mut r: Vec<usize> = Vec::new();
        for vi in &self.cells[&q].pts {
            r.push(*vi);
        }
        r
    }

    fn are_sibling_finalised(&self, qtc: &Vec<u8>) -> bool {
        let mut q2 = vec![0; qtc.len() - 1];
        q2.clone_from_slice(&qtc[..qtc.len() - 1]);
        q2.push(0);
        if self.cells.contains_key(&q2) && self.cells[&q2].finalised == true {
            q2.pop();
            q2.push(1);
            if self.cells.contains_key(&q2) && self.cells[&q2].finalised == true {
                q2.pop();
                q2.push(2);
                if self.cells.contains_key(&q2) && self.cells[&q2].finalised == true {
                    q2.pop();
                    q2.push(3);
                    if self.cells.contains_key(&q2) && self.cells[&q2].finalised == true {
                        return true;
                    }
                }
            }
        }
        false
    }

    fn set_qtc_as_finalised(&mut self, qtc: &Vec<u8>, b: bool) {
        self.cells.get_mut(qtc).unwrap().set_final(true);
    }

    fn get_cell_bbox(&self, gx: usize, gy: usize, gbbox: &mut [f64]) {
        gbbox[0] = self.minx + (gx * self.cellsize) as f64;
        gbbox[1] = self.miny + (gy * self.cellsize) as f64;
        gbbox[2] = self.minx + ((gx + 1) * self.cellsize) as f64;
        gbbox[3] = self.miny + ((gy + 1) * self.cellsize) as f64;
    }

    fn get_cell_bbox_qtc(&self, qtc: &Vec<u8>, gbbox: &mut [f64]) {
        let (mingx, mingy) = self.qtc2gxgy(qtc);
        let a: usize = (self.depth as usize) - qtc.len();
        let shift = 2_usize.pow(a as u32);
        gbbox[0] = self.minx + (mingx * self.cellsize) as f64;
        gbbox[1] = self.miny + (mingy * self.cellsize) as f64;
        gbbox[2] = self.minx + ((mingx + shift) * self.cellsize) as f64;
        gbbox[3] = self.miny + ((mingy + shift) * self.cellsize) as f64;
    }

    fn get_gxgy(&self, x: f64, y: f64) -> (usize, usize) {
        (
            ((x - self.minx) / self.cellsize as f64) as usize,
            ((y - self.miny) / self.cellsize as f64) as usize,
        )
    }

    fn get_qtc_from_gxgy(&self, gx: usize, gy: usize) -> Vec<u8> {
        let mut re: Vec<u8> = Vec::new();
        if self.depth == 0 {
            return re;
        }
        let mut mask: usize = 2_usize.pow(self.depth - 1);
        for i in 0..self.depth {
            let a = gx & mask == 2_usize.pow(self.depth - i - 1);
            let b = gy & mask == 2_usize.pow(self.depth - i - 1);
            if a == false && b == false {
                re.push(0)
            } else if a == false && b == true {
                re.push(1)
            } else if a == true && b == false {
                re.push(2)
            } else {
                re.push(3)
            }
            mask = mask >> 1;
        }
        re
    }

    fn get_qtc_from_xy(&self, x: f64, y: f64) -> Vec<u8> {
        let g = self.get_gxgy(x, y);
        self.get_qtc_from_gxgy(g.0, g.1)
    }

    fn qtc2gxgy(&self, qtc: &Vec<u8>) -> (usize, usize) {
        let mut q2 = vec![0; qtc.len()];
        q2.clone_from_slice(&qtc);
        self.gtc2gxgy_recursion(&qtc, 0, 0, 0)
    }

    fn gtc2gxgy_recursion(
        &self,
        c: &[u8],
        curdepth: usize,
        gx: usize,
        gy: usize,
    ) -> (usize, usize) {
        // let a: usize = (self.depth as usize) - c.len();
        let shift = self.griddim / (2_usize.pow(curdepth as u32)) / 2;
        if c.len() > 1 {
            if c[0] == 0 {
                return self.gtc2gxgy_recursion(&c[1..], curdepth + 1, gx, gy);
            } else if c[0] == 1 {
                return self.gtc2gxgy_recursion(&c[1..], curdepth + 1, gx, gy + shift);
            } else if c[0] == 2 {
                return self.gtc2gxgy_recursion(&c[1..], curdepth + 1, gx + shift, gy);
            } else {
                return self.gtc2gxgy_recursion(&c[1..], curdepth + 1, gx + shift, gy + shift);
            }
        }
        if c.is_empty() || c[0] == 0 {
            (gx, gy)
        } else if c[0] == 1 {
            (gx, gy + shift)
        } else if c[0] == 2 {
            (gx + shift, gy)
        } else {
            (gx + shift, gy + shift)
        }
    }
}

//----------------------
pub struct Triangulation {
    vs: Vec<[f64; 3]>,
    vs_no_incident: Vec<usize>, //-- number of incident triangles to each vertex (to know when to finalise)
    vs_active: Vec<bool>,       //-- is the vertex active (ie written once to the stream)
    vs_dead: Vec<bool>,         //-- is the vertex finalised (dead ids are swept from the qt cells at merges)
    vs_incident_tr: Vec<usize>, //-- link to one incident triangle
    ts: Vec<[usize; 6]>,
    qt: Quadtree,
    snaptol: f64,
    curt: usize,
    is_init: bool,
    robust_predicates: bool,
    freelist_ts: HashSet<usize>,
    freelist_vs: Vec<usize>,
    recent_ts: Vec<usize>,   //-- ring of the most recently created triangles
    recent_ts_cursor: usize, //-- position in the ring to write the next id
    sma_ids: HashMap<usize, usize>,
    sma_count_vertices: usize,
    out: io::BufWriter<io::Stdout>,
    pub max_vs: usize,
    pub max_ts: usize,
    pub n_inserts: u64,
    pub n_flips22: u64,
    pub n_walks: u64,
    pub n_walksteps: u64,
}

impl Triangulation {
    //-- new
    pub fn new() -> Triangulation {
        let mut thevs: Vec<[f64; 3]> = Vec::new();
        let p: [f64; 3] = [-99999.9, -99999.9, -99999.9];
        thevs.push(p); //-- add the infinity point
        let mut thevs_incident: Vec<usize> = Vec::new();
        thevs_incident.push(0);
        let mut thevs_i_tr: Vec<usize> = Vec::new();
        thevs_i_tr.push(0);
        let mut thevs_active: Vec<bool> = Vec::new();
        thevs_active.push(false);
        let mut thevs_dead: Vec<bool> = Vec::new();
        thevs_dead.push(true); //-- the infinity vertex is never a data point
        let mut thets: Vec<[usize; 6]> = Vec::new();
        let t: [usize; 6] = [0, 0, 0, 0, 0, 0]; //-- add a dummy triangle to simulate missing triangles (null pointers)
        thets.push(t);
        let q = Quadtree::new();
        let theflt: HashSet<usize> = HashSet::new();
        let theflv: Vec<usize> = Vec::new();
        let therecent: Vec<usize> = Vec::with_capacity(RECENT_TS_CAP);
        let thesma_ids: HashMap<usize, usize> = HashMap::new();
        Triangulation {
            vs: thevs,
            vs_no_incident: thevs_incident,
            vs_active: thevs_active,
            vs_dead: thevs_dead,
            vs_incident_tr: thevs_i_tr,
            ts: thets,
            qt: q,
            snaptol: 0.001,
            curt: 0,
            is_init: false,
            robust_predicates: true,
            freelist_vs: theflv,
            freelist_ts: theflt,
            recent_ts: therecent,
            recent_ts_cursor: 0,
            sma_ids: thesma_ids,
            sma_count_vertices: 1,
            //-- buffered output: one syscall per finalised cell instead of
            //-- one per line (flushes happen in finalise_qt_leaf/root)
            out: io::BufWriter::with_capacity(1 << 16, io::stdout()),
            max_vs: 0,
            max_ts: 0,
            n_inserts: 0,
            n_flips22: 0,
            n_walks: 0,
            n_walksteps: 0,
        }
    }

    pub fn set_cellsize(&mut self, c: usize) {
        self.qt.cellsize = c;
    }

    pub fn set_bbox(&mut self, minx: f64, miny: f64, maxx: f64, maxy: f64) {
        self.qt.minx = minx;
        self.qt.miny = miny;
        self.qt.maxx = maxx;
        self.qt.maxy = maxy;
    }

    pub fn set_grid_dimensions(&mut self, s: usize) {
        self.qt.init(s);
    }

    fn does_triangle_exist(&self, ti: usize, v0: usize, v1: usize, v2: usize) -> bool {
        let t0 = self.ts[ti];
        if t0[0] != v0 || t0[1] != v1 || t0[2] != v2 {
            false
        } else {
            true
        }
    }

    fn finalise_nonleaf_qtcell_and_merge(&mut self, qtc: &Vec<u8>) -> io::Result<()> {
        let mut q2 = vec![0; qtc.len() - 1];
        q2.clone_from_slice(&qtc[..qtc.len() - 1]);
        info!("Cell qtc{:?} finalised", q2);
        //-- get bbox of parent cell
        let mut gbbox: [f64; 4] = [0.0, 0.0, 0.0, 0.0];
        self.qt.get_cell_bbox_qtc(&q2, &mut gbbox);
        q2.push(0);
        //-- create parent cell
        let mut nc = Qtcell::new();
        nc.set_final(true);
        //-- copy from 4 children the vertices and triangles
        for i in 0..4 {
            q2.pop();
            q2.push(i);
            for (ti, vs) in &self.qt.cells[&q2].ts {
                if self.does_triangle_exist(*ti, vs.0, vs.1, vs.2) == true {
                    nc.add_ts(*ti, vs.0, vs.1, vs.2);
                }
            }
            for vi in &self.qt.cells[&q2].pts {
                if self.vs_dead[*vi] == false {
                    nc.add_pt(*vi);
                }
            }
            self.qt.cells.remove(&q2);
        }
        //-- insert the new parent cell in the QT
        q2.pop(); //-- now q2 is the new cell (the parent here)
        let mut req = vec![0; q2.len()];
        req.clone_from_slice(&q2);
        self.qt.cells.insert(req, nc);
        //-- verify if triangles encroach the new cell
        // let allts = self.qt.cells[&q2].get_vec_ts();
        let allts = &self.qt.cells.get(&q2).unwrap().ts;
        let mut fin_ts: Vec<usize> = Vec::new();
        for (ti, vs) in allts {
            if self.is_triangle_final(*ti, &gbbox) == true {
                fin_ts.push(*ti);
            }
        }
        for ti in &fin_ts {
            let _re = self.finalise_triangle(*ti, &q2);
        }
        //-- sweep the ids of dead vertices out of the new parent cell
        let vs_dead = &self.vs_dead;
        let q = self.qt.cells.get_mut(&q2).unwrap();
        q.pts.retain(|vi| vs_dead[*vi] == false);
        Ok(())
    }

    //-- finalise the leaf
    pub fn finalise_qt_leaf(&mut self, gx: usize, gy: usize) -> io::Result<()> {
        let mut qtc = self.qt.get_qtc_from_gxgy(gx, gy);
        info!(
            "Cell {}--{} finalised ({} vertices)",
            gx,
            gy,
            self.qt.cells[&qtc].number_pts()
        );
        let mut gbbox: [f64; 4] = [0.0, 0.0, 0.0, 0.0];
        self.qt.get_cell_bbox(gx, gy, &mut gbbox);

        //-- 1. collect all (unique) triangles and add to qt cell
        let allpts = &self.qt.cells.get(&qtc).unwrap().pts;
        let mut allts: HashMap<usize, (usize, usize, usize)> = HashMap::new();
        for vi in allpts {
            if self.vs_dead[*vi] == true {
                continue;
            }
            let l: Vec<usize> = self.incident_triangles_to_vertex(*vi).unwrap();
            for ti in &l {
                allts.insert(*ti, (self.ts[*ti][0], self.ts[*ti][1], self.ts[*ti][2]));
            }
        }

        //-- 2. check encroachment on current qt cell (a leaf at this point)
        //-- and delete them if final
        // let allts = &mut self.qt.cells.get_mut(&qtc).unwrap().ts;
        let mut t_to_keep: Vec<usize> = Vec::new();
        for (ti, vs) in &allts {
            // println!("{:?}", ti);
            if self.is_triangle_final(*ti, &gbbox) == true {
                let _re = self.finalise_triangle(*ti, &qtc);
            } else {
                t_to_keep.push(*ti);
            }
        }

        //-- 3. add the kept triangles to the qc cell
        let qc = &mut self.qt.cells.get_mut(&qtc).unwrap();
        for each in &t_to_keep {
            qc.add_ts(*each, allts[each].0, allts[each].1, allts[each].2);
        }

        //-- merge cell with parent, and continue until impossible
        //-- or only one cell left
        self.qt.set_qtc_as_finalised(&qtc, true);
        loop {
            if self.qt.are_sibling_finalised(&qtc) == true {
                let _re = self.finalise_nonleaf_qtcell_and_merge(&qtc);
                qtc.pop();
                if qtc.is_empty() {
                    let _re = self.finalise_qt_root();
                    break;
                }
            } else {
                break;
            }
        }
        //-- the output of this finalisation is complete: flush the
        //-- buffer so a downstream tool can consume it now
        self.out.flush()?;
        //-- sweep the ids of dead vertices out of this cell's points
        //-- (they were only marked as dead); doing this here —and not
        //-- only at the merge with the siblings— keeps the memory bound
        //-- tight even for incoherent streams where a cell waits for
        //-- its siblings until the very end
        let vs_dead = &self.vs_dead;
        let q = self.qt.cells.get_mut(&qtc).unwrap();
        q.pts.retain(|vi| vs_dead[*vi] == false);
        Ok(())
    }

    fn activate_vertex(&mut self, vi: usize) -> io::Result<()> {
        self.sma_ids.insert(vi, self.sma_count_vertices);
        self.sma_count_vertices += 1;
        write!(
            self.out,
            "v {} {} {}\n",
            &self.vs[vi][0], &self.vs[vi][1], &self.vs[vi][2]
        )?;
        self.vs_active[vi] = true;
        Ok(())
    }

    fn finalise_vertex(&mut self, vi: usize, qtc: &Vec<u8>) -> io::Result<()> {
        write!(self.out, "x {}\n", self.sma_ids[&vi])?;
        self.sma_ids.remove(&vi);
        //-- the vertex is only *marked* dead here; its id is swept from
        //-- the qt cell's list of points in one pass when the cell
        //-- merges with its siblings (or at the root). This keeps the
       //-- removal O(k) per cell instead of O(k²).
        self.vs_dead[vi] = true;
        self.freelist_vs.push(vi);
        Ok(())
    }

    /// write triangle to stdout stream + update neighbours + remove from qt cell
    fn finalise_triangle(&mut self, ti: usize, qtc: &Vec<u8>) -> io::Result<()> {
        let t0 = self.ts[ti];
        //-- check if the vertices involved are active, if not activate them
        for i in 0..3 {
            if self.vs_active[t0[i]] == false {
                let _re = self.activate_vertex(t0[i]);
            }
        }
        //-- write the triangle to the stream
        write!(
            self.out,
            "f {} {} {}\n",
            self.sma_ids[&t0[0]], self.sma_ids[&t0[1]], self.sma_ids[&t0[2]],
        )?;
        //-- add to free list to reuse in the future
        self.freelist_ts.insert(ti);
        // println!("{:?}", self.ts[ti]);
        //-- update neighbouring triangles (if any) for the walking
        for i in 3..6 {
            let tadj = t0[i];
            let k = &self.ts[tadj][3..6].iter().position(|&x| x == ti);
            if k.is_none() == false {
                self.ts[tadj][k.unwrap() + 3] = 0;
            }
        }
        //-- update # incident triangles for 3 vertices of the triangle
        //-- and finalise the vertex if now zero (never used again)
        for i in 0..3 {
            self.vs_no_incident[t0[i]] -= 1;
            if self.vs_no_incident[t0[i]] == 0 {
                let _re = self.finalise_vertex(t0[i], &qtc);
            }
        }
        //-- remove from qt cell list
        self.qt.cells.get_mut(qtc).unwrap().ts.remove(&ti);
        // println!("qtcell: {} -- {:?}", ti, self.qt.cells[qtc].ts);
        if self.curt == ti {
            //-- do not throw the walk hint away: re-point it to one of the
            //-- still-alive neighbours (eg. one straddling the frontier
            //-- between finalised and unfinalised space)
            self.curt = 0;
            for i in 3..6 {
                let n = t0[i];
                if n != 0
                    && self.freelist_ts.contains(&n) == false
                    && self.is_triangle_finite(n) == true
                {
                    self.curt = n;
                    break;
                }
            }
        }
        Ok(())
    }

    fn is_triangle_final(&self, ti: usize, bbox: &[f64]) -> bool {
        let t0 = self.ts[ti];
        //-- cannot remove triangles adjacent to boundary convexhull
        for i in 3..6 {
            if t0[i] != 0 && self.is_triangle_finite(t0[i]) == false {
                return false;
            }
        }
        //-- all 3 points inside the qt cell bbox?
        if geom::point_in_box(&self.vs[t0[0]], &bbox) == false
            || geom::point_in_box(&self.vs[t0[1]], &bbox) == false
            || geom::point_in_box(&self.vs[t0[2]], &bbox) == false
        {
            return false;
        }
        if geom::circumcircle_encroach_bbox(
            &self.vs[t0[0]],
            &self.vs[t0[1]],
            &self.vs[t0[2]],
            &bbox,
        ) == true
        {
            return false;
        }
        true
    }

    fn finalise_qt_root(&mut self) -> io::Result<()> {
        info!("Finalise the quadtree root cell");
        // let fout = format!("/Users/hugo/temp/sstout/leftovers.geojson");
        // let _re = self.write_geojson_triangles(fout.to_string());
        let q2: Vec<u8> = Vec::new();
        //-- triangles
        let allts = &self.qt.cells.get(&q2).unwrap().ts;
        let mut fin_ts: Vec<usize> = Vec::new();
        for (ti, vs) in allts {
            if self.is_triangle_finite(*ti) == true {
                fin_ts.push(*ti);
            }
        }
        for ti in &fin_ts {
            let _re = self.finalise_triangle(*ti, &q2);
        }

        let allpts = self.qt.cells.get(&q2).unwrap().pts.clone();
        for vi in &allpts {
            if self.vs_dead[*vi] == false {
                let _re = self.finalise_vertex(*vi, &q2);
            }
        }
        self.out.flush()?;
        Ok(())
    }

    //-- make sure everything written is pushed to the stream
    pub fn flush_output(&mut self) -> io::Result<()> {
        self.out.flush()
    }

    fn insert_one_pt_init_phase(&mut self, x: f64, y: f64, z: f64) -> Result<usize, usize> {
        let p: [f64; 3] = [x, y, z];
        for i in 0..self.vs.len() {
            if geom::distance2d_squared(&self.vs[i], &p) <= (self.snaptol * self.snaptol) {
                return Err(i);
            }
        }
        //-- add point to Triangulation
        let p = [x, y, z];
        // let c = self.vs.len();
        self.vs.push(p);
        self.vs_no_incident.push(0);
        self.vs_incident_tr.push(0);
        self.vs_active.push(false);
        self.vs_dead.push(false);
        let l = self.vs.len();
        if l >= 4 {
            let a = l - 3;
            let b = l - 2;
            let c = l - 1;
            let re = geom::orient2d(
                &self.vs[a],
                &self.vs[b],
                &self.vs[c],
                self.robust_predicates,
            );
            if re == 1 {
                // println!("init: ({},{},{})", a, b, c);

                self.ts.push([a, b, c, 3, 4, 2]);
                self.ts.push([a, 0, b, 3, 1, 4]);
                self.ts.push([b, 0, c, 4, 1, 2]);
                self.ts.push([c, 0, a, 2, 1, 3]);
                self.is_init = true;
            } else if re == -1 {
                // println!("init: ({},{},{})", a, c, b);
                self.ts.push([a, c, b, 3, 4, 2]);
                self.ts.push([a, 0, c, 3, 1, 4]);
                self.ts.push([c, 0, b, 4, 1, 2]);
                self.ts.push([b, 0, a, 2, 1, 3]);
                self.is_init = true;
            }
            self.vs_no_incident[a] = 3;
            self.vs_no_incident[b] = 3;
            self.vs_no_incident[c] = 3;
            //-- all 3 vertices point to the first real triangle ("1";
            //-- "0" is the dummy triangle, never a valid hint)
            self.vs_incident_tr[a] = 1;
            self.vs_incident_tr[b] = 1;
            self.vs_incident_tr[c] = 1;
        }
        self.curt = 1;
        // TODO: add those lines
        // if self.is_init == true {
        //     //-- insert the previous vertices in the dt
        //     for j in 1..(l - 3) {
        //         let tr = self.walk(&self.stars[&j].pt);
        //         // println!("found tr: {}", tr);
        //         self.flip13(j, &tr);
        //         self.update_dt(j);
        //     }
        // }
        Ok(self.vs.len() - 1)
    }

    /// Set a snap tolerance when inserting new points: if the newly inserted
    /// one is closer than snap_tolerance to another one, then it is not inserted.
    /// Avoids having very close vertices (like at 0.00007mm)
    /// Default is 0.001unit (thus 1mm for most datasets).
    pub fn set_snap_tolerance(&mut self, snaptol: f64) -> f64 {
        if snaptol > 0.0 {
            self.snaptol = snaptol;
        }
        self.snaptol
    }

    pub fn get_snap_tolerance(&self) -> f64 {
        self.snaptol
    }

    pub fn is_using_robust_predicates(&self) -> bool {
        self.robust_predicates
    }

    pub fn use_robust_predicates(&mut self, b: bool) {
        self.robust_predicates = b;
    }

    pub fn insert_one_pt_with_grid(&mut self, px: f64, py: f64, pz: f64) -> Result<usize, usize> {
        // println!("v: {}-{}", px, py);
        let re = self.insert_one_pt(px, py, pz);
        if re.is_ok() {
            let x = re.unwrap();
            self.qt.insert_one_vi(px, py, x);
        };
        //-- update for the stats
        if self.vs.len() > self.max_vs {
            self.max_vs = self.vs.len()
        }
        if self.ts.len() > self.max_ts {
            self.max_ts = self.ts.len()
        }
        re
    }

    //-- insert_one_pt
    fn insert_one_pt(&mut self, px: f64, py: f64, pz: f64) -> Result<usize, usize> {
        // println!("==>insertion ({}, {}, {})", px, py, pz);
        if self.is_init == false {
            return self.insert_one_pt_init_phase(px, py, pz);
        }
        //-- walk
        let p: [f64; 3] = [px, py, pz];
        let ti = self.walk(&p);
        // println!("Starting Ti: {}", ti);
        // assert_eq!(vec, [1, 2, 3]);
        //-- checking if within snap_tol of one vertex, then no insert
        for i in 0..3 {
            if geom::distance2d_squared(&self.vs[self.ts[ti][i]], &p)
                <= (self.snaptol * self.snaptol)
            {
                return Err(self.ts[ti][i]);
            }
        }
        //-- ok we now insert the point in the data structure
        let pi: usize;
        if self.freelist_vs.is_empty() {
            // println!("no free");
            pi = self.vs.len();
            self.vs.push([px, py, pz]);
            self.vs_no_incident.push(3); //-- b/c of flip13
            self.vs_active.push(false);
            self.vs_dead.push(false);
            self.vs_incident_tr.push(0);
        } else {
            // println!("freelist used: {:?}", self.freelist_vs);
            pi = self.freelist_vs.pop().unwrap();
            self.vs[pi] = [px, py, pz];
            self.vs_no_incident[pi] = 3; //-- b/c of flip13
            self.vs_active[pi] = false;
            self.vs_dead[pi] = false;
            self.vs_incident_tr[pi] = 0;
        }

        //-- flip13()
        let incidents = self.flip13(pi, ti);

        //-- add one of the new triangles as walk hint to the cell
        self.add_walk_hint_to_cell(pi, incidents);

        //-- update_dt()
        let mut mystack: Vec<usize> = Vec::new();
        mystack.push(incidents.0);
        mystack.push(incidents.1);
        mystack.push(incidents.2);

        loop {
            // self.print_ds();
            let ti = match mystack.pop() {
                None => break,
                Some(x) => x,
            };
            //-- fetch opposite
            let t = self.ts[ti];
            let k = &t[0..3].iter().position(|&x| x == pi).unwrap();
            let oppti = t[k + 3];
            if oppti == 0 {
                continue;
            }
            let oppt = self.ts[oppti];
            let ko = &oppt[3..6].iter().position(|&x| x == ti).unwrap();
            let oppvi = oppt[*ko];
            // assert_ne!(oppvi, 0);
            if oppvi == 0 {
                // println!("yo");
                continue;
            }
            if self.is_triangle_finite(ti) == false {
                let mut a: i8 = 0;
                if t[0] == 0 {
                    a = geom::orient2d(
                        &self.vs[oppvi],
                        &self.vs[t[1]],
                        &self.vs[t[2]],
                        self.robust_predicates,
                    );
                } else if t[1] == 0 {
                    a = geom::orient2d(
                        &self.vs[t[0]],
                        &self.vs[oppvi],
                        &self.vs[t[2]],
                        self.robust_predicates,
                    );
                } else if t[2] == 0 {
                    a = geom::orient2d(
                        &self.vs[t[0]],
                        &self.vs[t[1]],
                        &self.vs[oppvi],
                        self.robust_predicates,
                    );
                }
                if a > 0 {
                    self.flip22(ti, oppti);
                    mystack.push(ti);
                    mystack.push(oppti);
                }
            } else {
                if geom::incircle(
                    &self.vs[t[0]],
                    &self.vs[t[1]],
                    &self.vs[t[2]],
                    &self.vs[oppvi],
                    self.robust_predicates,
                ) > 0
                {
                    self.flip22(ti, oppti);
                    mystack.push(ti);
                    mystack.push(oppti);
                }
            }
        }
        //-- the self.curt is updated to one of the incident to the newly inserted point
        let mut newcurt = ti;
        if self.is_triangle_finite(newcurt) == false {
            let it = self.ts[newcurt];
            let k = &it[0..3].iter().position(|&x| x == 0).unwrap();
            newcurt = it[k + 3];
        }
        self.curt = newcurt;
        Ok(pi)
    }

    //-- remember a (re)created triangle in the "most recently created"
    //-- ring, used by the brute-force walk fallback
    fn add_recent_triangle(&mut self, ti: usize) {
        if self.recent_ts.len() < RECENT_TS_CAP {
            self.recent_ts.push(ti);
        } else {
            self.recent_ts[self.recent_ts_cursor] = ti;
        }
        self.recent_ts_cursor = (self.recent_ts_cursor + 1) % RECENT_TS_CAP;
    }

    fn flip13(&mut self, vi: usize, ti: usize) -> (usize, usize, usize) {
        self.n_inserts += 1;
        let newi0;
        let newi1;
        let mut reuse: bool = false;
        if self.freelist_ts.len() < 2 {
            newi0 = self.ts.len();
            newi1 = newi0 + 1;
        } else {
            // newi0 = self.freelist_ts.pop().unwrap();
            // newi1 = self.freelist_ts.pop().unwrap();
            newi0 = self.freelist_ts.iter().next().unwrap().clone();
            self.freelist_ts.remove(&newi0);
            newi1 = self.freelist_ts.iter().next().unwrap().clone();
            self.freelist_ts.remove(&newi1);
            reuse = true;
        }

        let t = self.ts[ti];
        let t0 = [vi, t[0], t[1], t[5], newi0, newi1];
        let t1 = [vi, t[1], t[2], t[3], newi1, ti];
        let t2 = [vi, t[2], t[0], t[4], ti, newi0];
        //-- increment the number of incident triangles
        self.vs_no_incident[t[0]] += 1;
        self.vs_no_incident[t[1]] += 1;
        self.vs_no_incident[t[2]] += 1;
        //-- set the pointer to one incident triangle to the new triangles
        self.vs_incident_tr[vi] = ti;
        self.vs_incident_tr[t[0]] = ti;
        self.vs_incident_tr[t[1]] = newi0;
        self.vs_incident_tr[t[2]] = newi1;

        //-- update neighbours
        let mut n = t[3];
        if n != 0 {
            let a = self.ts[n];
            let a1 = &a[3..6].iter().position(|&x| x == ti).unwrap();
            if let Some(x) = self.ts.get_mut(n) {
                (*x)[*a1 + 3] = newi0;
            }
        }
        n = t[4];
        if n != 0 {
            let a = self.ts[n];
            let a1 = &a[3..6].iter().position(|&x| x == ti).unwrap();
            if let Some(x) = self.ts.get_mut(n) {
                (*x)[*a1 + 3] = newi1;
            }
        }
        n = t[5];
        if n != 0 {
            let a = self.ts[n];
            let a1 = &a[3..6].iter().position(|&x| x == ti).unwrap();
            if let Some(x) = self.ts.get_mut(n) {
                (*x)[*a1 + 3] = ti;
            }
        }
        if let Some(x) = self.ts.get_mut(ti) {
            (*x) = t0;
        }
        if reuse == false {
            self.ts.push(t1);
            self.ts.push(t2);
        } else {
            self.ts[newi0] = t1;
            self.ts[newi1] = t2;
        }
        //-- remember the 3 new triangles for the brute-force walk fallback
        self.add_recent_triangle(ti);
        self.add_recent_triangle(newi0);
        self.add_recent_triangle(newi1);
        (ti, newi0, newi1)
    }

    pub fn incident_triangles_to_vertex(&self, vi: usize) -> Option<Vec<usize>> {
        let mut l: Vec<usize> = Vec::new();
        // let t0 = self.walk(&[self.vs[vi][0], self.vs[vi][1]]);
        let t0 = self.vs_incident_tr[vi];
        // println!("{}--{:?}", t0, self.ts[t0]);
        l.push(t0);
        let mut curt = t0;
        loop {
            let tt = self.ts[curt];
            let k0 = &tt[0..3].iter().position(|&x| x == vi).unwrap();
            let nextt = tt[((*k0 + 1) % 3) + 3];
            // println!("cur: {}", nextt);
            if nextt != t0 {
                l.push(nextt);
                curt = nextt;
            } else {
                break;
            }
        }
        // println!("{:?}", l);
        Some(l)
    }

    fn is_triangle_finite(&self, t: usize) -> bool {
        let t = self.ts[t];
        if t[0] == 0 || t[1] == 0 || t[2] == 0 {
            false
        } else {
            true
        }
    }

    /// Returns number of finite vertices in the triangulation.
    pub fn number_of_vertices(&self) -> usize {
        //-- number of finite vertices
        self.vs.len() - self.freelist_vs.len() - 1
    }

    /// Returns number of finite triangles in the triangulation.
    fn number_of_triangles(&self) -> usize {
        //-- number of finite vertices
        let mut total: usize = 0;
        for (i, t) in self.ts.iter().enumerate() {
            if self.is_triangle_finite(i) == true {
                total += 1;
            }
        }
        total
    }

    //-- add one of the triangles incident to the newly inserted point as
    //-- a restart hint to the quadtree cell containing that point
    //-- (mechanism of Isenburg et al. §4.1: a walk that falls into
    //-- finalised space is restarted from triangles known to be near)
    fn add_walk_hint_to_cell(&mut self, vi: usize, tris: (usize, usize, usize)) {
        if vi == 0 {
            return;
        }
        //-- prefer a finite triangle as hint (an infinite one would be
        //-- rejected by the walk anyway)
        let mut thet: usize = tris.0;
        if self.is_triangle_finite(thet) == false {
            thet = tris.1;
            if self.is_triangle_finite(thet) == false {
                thet = tris.2;
            }
        }
        let g = self.qt.get_gxgy(self.vs[vi][0], self.vs[vi][1]);
        let qtc = self.qt.get_qtc_from_gxgy(g.0, g.1);
        if let Some(c) = self.qt.cells.get_mut(&qtc) {
            if c.finalised == false {
                c.add_walk_hint(thet);
            }
        }
    }

    fn walk(&mut self, x: &[f64]) -> usize {
        self.n_walks += 1;
        //-- 1. try walk from latest
        let tr = self.curt;
        if tr != 0 {
            let re = self.walk_safe(x, tr);
            if re.is_some() {
                return re.unwrap();
            }
        }

        //-- 2. try to restart the walk from the triangles hinted in the
        //--    cell containing the point, then in its neighbouring cells
        //--    (cf. Isenburg et al. §4.1: restart from triangles on the
        //--    quadrant's list, then on neighbouring quadrants' lists)
        let g = self.qt.get_gxgy(x[0], x[1]);
        let qtc = self.qt.get_qtc_from_gxgy(g.0, g.1);
        let mut cells_to_try: Vec<Vec<u8>> = Vec::new();
        cells_to_try.push(qtc);
        for dy in -1..2 as isize {
            for dx in -1..2 as isize {
                if dx == 0 && dy == 0 {
                    continue; //-- own cell is already first in the list
                }
                let gx = g.0 as isize + dx;
                let gy = g.1 as isize + dy;
                if gx < 0 || gy < 0 {
                    continue;
                }
                let gx = gx as usize;
                let gy = gy as usize;
                if gx >= self.qt.griddim || gy >= self.qt.griddim {
                    continue;
                }
                let qtc2 = self.qt.get_qtc_from_gxgy(gx, gy);
                cells_to_try.push(qtc2);
            }
        }
        let mut hint_candidates: Vec<usize> = Vec::new();
        for eachqtc in &cells_to_try {
            if let Some(cell) = self.qt.cells.get(eachqtc) {
                for each in &cell.walk_hints {
                    if *each != 0
                        && self.freelist_ts.contains(each) == false
                        && self.is_triangle_finite(*each) == true
                    {
                        hint_candidates.push(*each);
                    }
                }
            }
        }
        for each in &hint_candidates {
            let re = self.walk_safe(x, *each);
            if re.is_some() {
                return re.unwrap();
            }
        }

        //-- 3. try from a few randomly chosen vertices of the cell
        //--    (bounded: at most 256 attempts, each with a valid hint)
        let empty: Vec<usize> = Vec::new();
        let allpts: Vec<usize> = match self.qt.cells.get(&cells_to_try[0]) {
            Some(c) => c.pts.clone(),
            None => empty,
        };
        if allpts.is_empty() == false {
            warn!("walk.3");
            let mut rng = thread_rng();
            let nattempts = allpts.len().min(256);
            let mut candidates: Vec<usize> = Vec::with_capacity(nattempts);
            for _i in 0..nattempts {
                let vi = allpts[rng.gen_range(0..allpts.len())];
                if self.vs_dead[vi] == true {
                    continue;
                }
                let t0 = self.vs_incident_tr[vi];
                if t0 != 0
                    && self.freelist_ts.contains(&t0) == false
                    && self.is_triangle_finite(t0) == true
                {
                    candidates.push(t0);
                }
            }
            for t0 in &candidates {
                let re = self.walk_safe(x, *t0);
                if re.is_some() {
                    return re.unwrap();
                }
            }
        }

        //-- 3.5. the walk cannot reach the containing triangle (there is
        //--      finalised space in between), but for a query in an
        //--      unfinalised cell it is incident to one of that cell's
        //--      vertices: search the stars of the cell's points directly
        if allpts.is_empty() == false {
            for vi in &allpts {
                if self.vs_dead[*vi] == true {
                    continue;
                }
                let l = self.incident_triangles_to_vertex(*vi).unwrap();
                for ti in &l {
                    if self.freelist_ts.contains(ti) == false
                        && self.is_triangle_finite(*ti) == true
                    {
                        let t = self.ts[*ti];
                        if geom::intriangle(
                            &self.vs[t[0]],
                            &self.vs[t[1]],
                            &self.vs[t[2]],
                            &x,
                            self.robust_predicates,
                        ) == 1
                        {
                            return *ti;
                        }
                    }
                }
            }
        }

        //-- 4. try brute-force
        //-- TODO: this brute-force too?
        // let re2 = self.walk_bruteforce_closest_vertex_then_walksafe(x);
        // if re2.is_some() {
        //     return re2.unwrap();
        // }
        warn!("walk.4");
        let re3 = self.walk_bruteforce_triangles(x);
        if re3.is_some() {
            return re3.unwrap();
        }

        //-- 4. we are outside the CH of the current dataset
        // warn!("point is outside the CH, finding closest point on the CH");
        warn!("walk.5");
        let re4 = self.walk_bruteforce_outside_convex_hull(x);
        if re4.is_some() {
            return re4.unwrap();
        } else {
            error!("WALK FAILED MISERABLY :'(");
        }

        return 0;
    }

    fn walk_safe(&mut self, x: &[f64], mut tr: usize) -> Option<usize> {
        if tr == 0 {
            //-- the dummy triangle: never a valid starting point
            return None;
        }
        loop {
            self.n_walksteps += 1;
            if self.is_triangle_finite(tr) == false {
                return Some(tr);
            }
            let t = self.ts[tr];
            if geom::orient2d(&self.vs[t[0]], &self.vs[t[1]], &x, self.robust_predicates) != -1 {
                if geom::orient2d(&self.vs[t[1]], &self.vs[t[2]], &x, self.robust_predicates) != -1
                {
                    if geom::orient2d(&self.vs[t[2]], &self.vs[t[0]], &x, self.robust_predicates)
                        != -1
                    {
                        // println!("found it! {}", tr);
                        return Some(tr);
                    } else {
                        if t[4] == 0 {
                            //-- if about to fall off then return None
                            return None;
                        }
                        tr = t[4];
                    }
                } else {
                    if t[3] == 0 {
                        return None;
                    }
                    tr = t[3];
                }
            } else {
                if t[5] == 0 {
                    return None;
                }
                tr = t[5];
            }
        }
    }

    fn is_vertex_convex_hull(&self, vi: usize) -> bool {
        true
    }

    fn walk_bruteforce_outside_convex_hull(&self, x: &[f64]) -> Option<usize> {
        // warn!("walk brute-force ON CONVEX HULL");
        let mut amin: f64 = std::f64::MAX;
        let mut tmin: usize = 0;
        for (id, t) in self.ts.iter().enumerate() {
            if self.is_triangle_finite(id) == false {
                let t = self.ts[id];
                let a: f64;
                if t[0] == 0 {
                    a = geom::area_triangle(&x, &self.vs[t[1]], &self.vs[t[2]]);
                } else if t[1] == 0 {
                    a = geom::area_triangle(&self.vs[t[0]], &x, &self.vs[t[2]]);
                } else {
                    a = geom::area_triangle(&self.vs[t[0]], &self.vs[t[1]], &x);
                }
                if a > 0.0 && a < amin {
                    amin = a;
                    tmin = id;
                }
            }
        }
        Some(tmin)
    }

    fn walk_bruteforce_triangles(&self, x: &[f64]) -> Option<usize> {
        // warn!("walk_bruteforce_triangles()");
        //-- the triangle containing x is almost always one of the most
        //-- recently created ones: scan the ring first (newest first)
        let nring = self.recent_ts.len();
        let mut idx = self.recent_ts_cursor;
        for _i in 0..nring {
            idx = (idx + RECENT_TS_CAP - 1) % RECENT_TS_CAP;
            let ti = self.recent_ts[idx];
            if self.freelist_ts.contains(&ti) == false && self.is_triangle_finite(ti) == true {
                let t = self.ts[ti];
                if geom::intriangle(
                    &self.vs[t[0]],
                    &self.vs[t[1]],
                    &self.vs[t[2]],
                    &x,
                    self.robust_predicates,
                ) == 1
                {
                    return Some(ti);
                }
            }
        }
        //-- not found among the recent ones: full sweep
        for (id, t) in self.ts.iter().enumerate() {
            if self.freelist_ts.contains(&id) == false && self.is_triangle_finite(id) == true {
                if geom::intriangle(
                    &self.vs[t[0]],
                    &self.vs[t[1]],
                    &self.vs[t[2]],
                    &x,
                    self.robust_predicates,
                ) == 1
                {
                    return Some(id);
                }
            }
        }
        return None;
    }

    // fn walk_bruteforce_closest_vertex_then_walksafe(&self, x: &[f64]) -> Option<usize> {
    //     //-- find closest vertex that is on the CH
    //     let mut dmin: f64 = std::f64::MAX;
    //     let mut vmin: usize = 0;
    //     for i in self.stars.keys() {
    //         if *i != 0 {
    //             let d = geom::distance2d_squared(x, &self.vs[i]);
    //             if d < dmin {
    //                 dmin = d;
    //                 vmin = *i;
    //             }
    //         }
    //     }
    //     self.walk_safe(x, vmin)
    // }

    fn flip22(&mut self, t0i: usize, t1i: usize) {
        self.n_flips22 += 1;
        // println!("flip22");
        // self.print_ds();
        let t0 = self.ts[t0i];
        let t1 = self.ts[t1i];
        let k0 = &t0[3..6].iter().position(|&x| x == t1i).unwrap();
        let k1 = &t1[3..6].iter().position(|&x| x == t0i).unwrap();

        let nt0 = [
            t0[*k0],
            t0[(*k0 + 1) % 3],
            t1[*k1],
            t1[((*k1 + 1) % 3) + 3],
            t1i,
            t0[((*k0 + 2) % 3) + 3],
        ];
        let nt1 = [
            t0[*k0],
            t1[*k1],
            t0[(*k0 + 2) % 3],
            t1[((*k1 + 2) % 3) + 3],
            t0[((*k0 + 1) % 3) + 3],
            t0i,
        ];
        //-- update the number of incident triangles
        self.vs_no_incident[t0[*k0]] += 1;
        self.vs_no_incident[t1[*k1]] += 1;
        self.vs_no_incident[t0[(*k0 + 1) % 3]] -= 1;
        self.vs_no_incident[t0[(*k0 + 2) % 3]] -= 1;
        //-- update the pointer incident triangle for each 4 vertices involved
        self.vs_incident_tr[t0[*k0]] = t0i;
        self.vs_incident_tr[t1[*k1]] = t1i;
        self.vs_incident_tr[t0[(*k0 + 1) % 3]] = t0i;
        self.vs_incident_tr[t0[(*k0 + 2) % 3]] = t1i;

        //-- update 2 neighbouring triangles
        //-- because 2 stay the same since the indices are not changed
        let mut n0 = t0[((*k0 + 1) % 3) + 3];
        if n0 != 0 {
            let a = self.ts[n0];
            let a1 = &a[3..6].iter().position(|&x| x == t0i).unwrap();
            if let Some(x) = self.ts.get_mut(n0) {
                (*x)[*a1 + 3] = t1i;
            }
        }
        n0 = t1[((k1 + 1) % 3) + 3];
        if n0 != 0 {
            let a = self.ts[n0];
            let a1 = &a[3..6].iter().position(|&x| x == t1i).unwrap();
            if let Some(x) = self.ts.get_mut(n0) {
                (*x)[*a1 + 3] = t0i;
            }
        }
        if let Some(x) = self.ts.get_mut(t0i) {
            (*x) = nt0;
        }
        if let Some(x) = self.ts.get_mut(t1i) {
            (*x) = nt1;
        }
        //-- remember the 2 recreated triangles for the brute-force walk fallback
        self.add_recent_triangle(t0i);
        self.add_recent_triangle(t1i);
    }

    pub fn print_ds(&self) {
        println!("=== vertices ===");
        // let mut ids = self.vs.keys().copied().collect::<Vec<_>>();
        // ids.sort();
        for (id, v) in self.vs.iter().enumerate() {
            println!(
                "{} -- ({}, {}, {})",
                id, self.vs[id][0], self.vs[id][1], self.vs[id][2]
            );
        }
        println!("=== triangles ===");
        // ids = self.ts.keys().copied().collect::<Vec<_>>();
        // ids.sort();
        for (id, t) in self.ts.iter().enumerate() {
            println!(
                "{} -- [{}, {}, {}, {}, {}, {}]",
                id,
                self.ts[id][0],
                self.ts[id][1],
                self.ts[id][2],
                self.ts[id][3],
                self.ts[id][4],
                self.ts[id][5]
            );
        }
    }

    /// write a GeoJSON file of the quadtree/grid to disk
    pub fn write_geojson_grid(&self, path: String) -> std::io::Result<()> {
        let mut fc = FeatureCollection {
            bbox: None,
            features: vec![],
            foreign_members: None,
        };
        for (key, cell) in &self.qt.cells {
            // println!("key:{:?}", key);
            let mut bbox: [f64; 4] = [0.0, 0.0, 0.0, 0.0];
            self.qt.get_cell_bbox_qtc(key, &mut bbox);
            // println!("bbox: {:?}", bbox);
            let mut l: Vec<Vec<Vec<f64>>> = vec![vec![Vec::with_capacity(1); 5]];
            l[0][0].push(bbox[0]);
            l[0][0].push(bbox[1]);
            l[0][1].push(bbox[2]);
            l[0][1].push(bbox[1]);
            l[0][2].push(bbox[2]);
            l[0][2].push(bbox[3]);
            l[0][3].push(bbox[0]);
            l[0][3].push(bbox[3]);
            l[0][4].push(bbox[0]);
            l[0][4].push(bbox[1]);
            let mut attributes = Map::new();
            let n = self.qt.cells[key].number_pts();
            attributes.insert(String::from("no_pts"), to_value(n).unwrap());
            attributes.insert(String::from("qtc"), to_value(key).unwrap());
            let gtr = Geometry::new(Value::Polygon(l));
            let f = Feature {
                bbox: None,
                geometry: Some(gtr),
                id: None,
                properties: Some(attributes),
                foreign_members: None,
            };
            fc.features.push(f);
        }
        //-- write the file to disk
        let mut fo = File::create(path)?;
        write!(fo, "{}", fc.to_string()).unwrap();
        Ok(())
    }

    /// write a GeoJSON file of the triangles/vertices to disk
    pub fn write_geojson_triangles(&self, path: String) -> std::io::Result<()> {
        let mut fc = FeatureCollection {
            bbox: None,
            features: vec![],
            foreign_members: None,
        };
        //-- vertices
        for (i, v) in self.vs.iter().enumerate() {
            if i == 0 {
                continue;
            }
            let pt = Geometry::new(Value::Point(vec![v[0], v[1]]));
            let mut attributes = Map::new();
            attributes.insert(String::from("id"), to_value(i.to_string()).unwrap());
            attributes.insert(
                String::from("incident"),
                to_value(self.vs_no_incident[i].to_string()).unwrap(),
            );
            let f = Feature {
                bbox: None,
                geometry: Some(pt),
                id: None,
                properties: Some(attributes),
                foreign_members: None,
            };
            fc.features.push(f);
        }
        //-- triangles
        for (i, t) in self.ts.iter().enumerate() {
            if self.is_triangle_finite(i) == false {
                continue;
            }
            if self.freelist_ts.contains(&i) {
                continue;
            }
            let mut l: Vec<Vec<Vec<f64>>> = vec![vec![Vec::with_capacity(1); 4]];
            l[0][0].push(self.vs[t[0]][0]);
            l[0][0].push(self.vs[t[0]][1]);
            l[0][1].push(self.vs[t[1]][0]);
            l[0][1].push(self.vs[t[1]][1]);
            l[0][2].push(self.vs[t[2]][0]);
            l[0][2].push(self.vs[t[2]][1]);
            l[0][3].push(self.vs[t[0]][0]);
            l[0][3].push(self.vs[t[0]][1]);
            let gtr = Geometry::new(Value::Polygon(l));
            let mut attributes = Map::new();
            attributes.insert(String::from("tid"), to_value(i.to_string()).unwrap());
            let f = Feature {
                bbox: None,
                geometry: Some(gtr),
                id: None,
                properties: Some(attributes),
                foreign_members: None,
            };
            fc.features.push(f);
        }
        //-- write the file to disk
        let mut fo = File::create(path)?;
        write!(fo, "{}", fc.to_string()).unwrap();
        Ok(())
    }
}

impl fmt::Display for Triangulation {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        fmt.write_str("======== TRIANGULATION ========\n")?;
        fmt.write_str(&format!("# vertices: {:19}\n", self.number_of_vertices()))?;
        fmt.write_str(&format!("# triangles: {:18}\n", self.number_of_triangles()))?;
        // fmt.write_str(&format!(
        //     "# convex hull: {:16}\n",
        //     self.number_of_vertices_on_convex_hull()
        // ))?;
        fmt.write_str(&format!("---\nrobust: {}\n", self.robust_predicates))?;
        fmt.write_str("===============================\n")?;
        Ok(())
    }
}
