//! # startin

pub mod geom;

use std::fmt;
use std::fs::File;
use std::io::Write;
use std::io::{self};

use geojson::{Feature, FeatureCollection, Geometry, Value};
use serde_json::{to_value, Map};

use hashbrown::HashMap;
use std::collections::HashSet;

extern crate rand;

//----------------------
struct Qtcell {
    pts: HashSet<usize>,
    ts: HashSet<usize>,
    finalised: bool,
}

impl Qtcell {
    pub fn new() -> Qtcell {
        let p: HashSet<usize> = HashSet::new();
        let t: HashSet<usize> = HashSet::new();
        Qtcell {
            pts: p,
            ts: t,
            finalised: false,
        }
    }

    pub fn number_pts(&self) -> usize {
        self.pts.len()
    }

    pub fn number_ts(&self) -> usize {
        self.ts.len()
    }

    pub fn add_pt(&mut self, vi: usize) {
        self.pts.insert(vi);
    }

    pub fn add_ts(&mut self, ti: usize) {
        self.ts.insert(ti);
    }

    pub fn set_final(&mut self, b: bool) {
        self.finalised = b
    }

    pub fn get_vec_vs(&self) -> Vec<usize> {
        let mut r: Vec<usize> = Vec::new();
        for vi in &self.pts {
            r.push(*vi);
        }
        r
    }
    pub fn get_vec_ts(&self) -> Vec<usize> {
        let mut r: Vec<usize> = Vec::new();
        for ti in &self.ts {
            r.push(*ti);
        }
        r
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

    pub fn get_cell_count(&self, gx: usize, gy: usize) -> Option<usize> {
        if gx >= self.griddim || gy > self.griddim {
            None
        } else {
            let q = self.get_qtc_from_gxgy(gx, gy);
            Some(self.cells[&q].pts.len())
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

    fn are_sibling_final(&self, qtc: &Vec<u8>) -> bool {
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

    fn set_final_qtc(&mut self, qtc: &Vec<u8>, b: bool) {
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
    ts: Vec<[usize; 6]>,
    qt: Quadtree,
    snaptol: f64,
    curt: usize,
    is_init: bool,
    robust_predicates: bool,
    freelist_ts: Vec<usize>, // TODO: use the freelist for insertion
    freelist_vs: Vec<usize>,
    sma_ids: HashMap<usize, usize>,
    sma_count_vertices: usize,
}

impl Triangulation {
    //-- new
    pub fn new() -> Triangulation {
        let mut thevs: Vec<[f64; 3]> = Vec::new();
        let p: [f64; 3] = [-99999.9, -99999.9, -99999.9];
        thevs.push(p); //-- add the infinity point
        let mut thets: Vec<[usize; 6]> = Vec::new();
        let t: [usize; 6] = [0, 0, 0, 0, 0, 0]; //-- add a dummy triangle to simulate missing triangles (null pointers)
        thets.push(t);
        let q = Quadtree::new();
        let theflt: Vec<usize> = Vec::new();
        let theflv: Vec<usize> = Vec::new();
        let thesma_ids: HashMap<usize, usize> = HashMap::new();
        Triangulation {
            vs: thevs,
            ts: thets,
            qt: q,
            snaptol: 0.001,
            curt: 0,
            is_init: false,
            robust_predicates: true,
            freelist_vs: theflv,
            freelist_ts: theflt,
            sma_ids: thesma_ids,
            sma_count_vertices: 1,
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

    fn finalise_cell_and_merge(&mut self, qtc: &Vec<u8>) -> io::Result<()> {
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
        //-- copy from 4 children the vertices/triangles that remain active
        let mut active_vs_in_ts: HashSet<usize> = HashSet::new();
        let mut all_vs: HashSet<usize> = HashSet::new();
        for i in 0..4 {
            q2.pop();
            q2.push(i);
            let allts: Vec<usize> = self.qt.cells[&q2].get_vec_ts();
            for ti in &allts {
                if self.is_triangle_final(*ti, &gbbox) == true {
                    let _re = self.finalise_triangle(*ti);
                } else {
                    nc.add_ts(*ti);
                    for i in 0..3 {
                        active_vs_in_ts.insert(self.ts[*ti][i]);
                    }
                }
            }
            for vi in &self.qt.cells[&q2].pts {
                all_vs.insert(*vi);
            }
            self.qt.cells.remove(&q2);
        }
        //-- keep only the active vertices and add them to nc
        info!("all_vs: {:?}", all_vs);
        info!("active_vs_in_ts: {:?}", active_vs_in_ts);
        let intersect: HashSet<_> = all_vs.intersection(&active_vs_in_ts).collect();
        info!("still active (intersection): {:?}", intersect);
        for vi in &intersect {
            nc.add_pt(**vi);
        }
        let difference: HashSet<_> = all_vs.difference(&active_vs_in_ts).collect();
        info!("finalised (difference): {:?}", difference);
        for vi in &difference {
            io::stdout().write_all(&format!("x {}\n", self.sma_ids[*vi]).as_bytes())?;
        }
        //-- insert the new parent cell in the QT
        q2.pop();
        let mut req = vec![0; q2.len()];
        req.clone_from_slice(&q2);
        self.qt.cells.insert(q2, nc);
        Ok(())
    }

    pub fn finalise_qtcell(&mut self, gx: usize, gy: usize) -> io::Result<()> {
        info!(
            "Cell {}--{} finalised ({} vertices)",
            gx,
            gy,
            self.qt.get_cell_count(gx, gy).unwrap()
        );

        let mut qtc = self.qt.get_qtc_from_gxgy(gx, gy);

        let mut gbbox: [f64; 4] = [0.0, 0.0, 0.0, 0.0];
        self.qt.get_cell_bbox(gx, gy, &mut gbbox);

        //-- 1. collect all (unique) triangles
        // let mut allts: HashSet<usize> = HashSet::new();
        let mut allts: HashMap<usize, bool> = HashMap::new();
        let allpts: Vec<usize> = self.qt.get_cell_pts(gx, gy);

        //-- active all the vertices of the leaf cell
        //-- TODO: maybe make this better (benchmark?)
        for (i, vi) in allpts.iter().enumerate() {
            self.sma_ids.insert(*vi, self.sma_count_vertices);
            self.sma_count_vertices += 1;
            io::stdout().write_all(
                &format!(
                    "v {} {} {}\n",
                    &self.vs[*vi][0], &self.vs[*vi][1], &self.vs[*vi][2]
                )
                .as_bytes(),
            )?;
        }
        //-- check the star of each vertex and collect triangles
        //-- and finalise the vertices whose star is filled with triangles
        //-- that can be finalised
        let mut fin_vs: Vec<usize> = Vec::new();
        for vi in &allpts {
            //-- 2. encroachment?
            let l: Vec<usize> = self.incident_triangles_to_vertex(*vi).unwrap();
            let mut finv: bool = true;
            for ti in &l {
                if allts.contains_key(ti) {
                    if allts[ti] == false {
                        finv = false;
                    }
                } else if self.is_triangle_final(*ti, &gbbox) == false {
                    finv = false;
                    allts.insert(*ti, false);
                } else {
                    // //-- write to stream final triangles
                    // io::stdout().write_all(
                    //     &format!(
                    //         "f {} {} {}\n",
                    //         self.sma_ids[&self.ts[*ti][0]],
                    //         self.sma_ids[&self.ts[*ti][1]],
                    //         self.sma_ids[&self.ts[*ti][2]],
                    //     )
                    //     .as_bytes(),
                    // )?;
                    allts.insert(*ti, true);
                }
            }
            if finv == true {
                fin_vs.push(*vi);
                // println!("vertex final {:?}", *vi);
                //-- TODO: remove the point from the DS!
                //-- finalise the vertices in the stream
                // io::stdout().write_all(&format!("x {}\n", self.sma_ids[vi]).as_bytes())?;
                // self.qt.cells.get_mut(&qtc).unwrap().pts.remove(vi);
            }
        }
        // println!("allts: {:?}", allts);

        //-- finalise the triangles (write to stream)
        //-- and remove from the ds
        for (ti, finalised) in &allts {
            // println!("{}: \"{}\"", ti, finalised);
            if *finalised == true {
                let _re = self.finalise_triangle(*ti);
            } else {
                self.qt.cells.get_mut(&qtc).unwrap().add_ts(*ti);
            }
        }

        //-- finalise the vertices that will never be used again
        for vi in &fin_vs {
            //-- TODO: remove the point from the DS!
            //-- finalise the vertices in the stream
            io::stdout().write_all(&format!("x {}\n", self.sma_ids[vi]).as_bytes())?;
            self.qt.cells.get_mut(&qtc).unwrap().pts.remove(vi);
        }

        //-- merge cell with parent, and continue until impossible
        //-- or only one cell left
        self.qt.set_final_qtc(&qtc, true);
        loop {
            if self.qt.are_sibling_final(&qtc) == true {
                self.finalise_cell_and_merge(&qtc);
                qtc.pop();
                if qtc.is_empty() {
                    self.finalise_qt_root();
                    break;
                }
            } else {
                break;
            }
        }
        Ok(())
    }

    fn finalise_triangle(&mut self, ti: usize) -> io::Result<()> {
        io::stdout().write_all(
            &format!(
                "f {} {} {}\n",
                self.sma_ids[&self.ts[ti][0]],
                self.sma_ids[&self.ts[ti][1]],
                self.sma_ids[&self.ts[ti][2]],
            )
            .as_bytes(),
        )?;
        self.freelist_ts.push(ti);
        // println!("{:?}", self.ts[ti]);
        for i in 3..6 {
            let tadj = self.ts[ti][i];
            let k = &self.ts[tadj][3..6].iter().position(|&x| x == ti);
            if k.is_none() == false {
                self.ts[tadj][k.unwrap() + 3] = 0;
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

    // let qtc = self.qt.finalise_cell(gx, gy);

    // if (qtc.is_empty() == true) || (self.number_of_vertices() <= 3) {
    //     return Ok(());
    // }
    // //-- nothing to do if single cell finalised and it's empty
    // if (qtc.len() as u32 == self.qt.depth) && (self.qt.get_cell_count(gx, gy).unwrap() == 0) {
    //     return Ok(());
    // }

    // let allfcells = self.qt.get_all_gxgy_from_qtc(&qtc);

    // let mut gbbox: [f64; 4] = [0.0, 0.0, 0.0, 0.0];
    // self.qt.get_cell_bbox_qtc(&qtc, &mut gbbox);
    // for c in allfcells {
    //     let allvinc = self.qt.gpts[c.0][c.1].clone();
    //     for theid in allvinc {
    //         let mut re = self.adjacent_vertices_to_vertex(theid).unwrap();
    //         let mut repts: Vec<Option<Vec<f64>>> = Vec::with_capacity(re.len());

    //         for each in &re {
    //             repts.push(self.get_point(*each));
    //         }
    //         //-- check if edge on CH
    //         let mut onch = 0;
    //         for (i, _p) in re.iter().enumerate() {
    //             if repts[i].is_some() && self.stars[&re[i]].link.contains_infinite_vertex() {
    //                 onch += 1;
    //             }
    //         }
    //         if onch >= 2 {
    //             continue;
    //         }
    //         let theidpt = self.get_point(theid).unwrap();
    //         let mut fin: bool = true;
    //         re.push(re[0]);
    //         repts.push(repts[0].clone());
    //         for i in 0..(re.len() - 1) {
    //             if repts[i].is_some() == true
    //                 && repts[i + 1].is_some() == true
    //                 && geom::circumcentre_encroach_bbox(
    //                     &theidpt,
    //                     &repts[i].as_ref().unwrap(),
    //                     &repts[i + 1].as_ref().unwrap(),
    //                     &gbbox,
    //                 ) == true
    //             {
    //                 fin = false;
    //                 break;
    //             }
    //         }
    //         re.pop();
    //         if fin == true {
    //             match self.outputmode {
    //                 // write triangles
    //                 Outputmode::Sma => {
    //                     self.write_sma_one_vertex(theid)
    //                         .expect("Failed to write vertex");
    //                     //-- flush it from QT and the DS
    //                     self.qt.gpts[c.0][c.1].remove(&theid);
    //                     self.flush_star(theid);
    //                 }

    //                 //-- finalise the vertex with its star/link
    //                 Outputmode::Stars => {
    //                     io::stdout()
    //                         .write_all(&format!("x {} {:?}\n", theid, re).as_bytes())?;
    //                     //-- flush it from QT and the DS
    //                     self.qt.gpts[c.0][c.1].remove(&theid);
    //                     self.flush_star(theid);
    //                 }

    //                 //-- finalise the vertex, output both vertex info and star data
    //                 Outputmode::Both => {
    //                     self.write_stars_one_vertex(theid)
    //                         .expect("Failed to write vertex");
    //                     //-- flush it from QT and the DS
    //                     self.qt.gpts[c.0][c.1].remove(&theid);
    //                     self.flush_star(theid);
    //                 }

    //                 _ => {
    //                     println!("Option not implemented");
    //                 }
    //             }
    //         }
    //     }
    // }

    // fn write_sma_one_vertex(&mut self, v: usize) -> io::Result<()> {
    //     let vpt = self.get_point(v).unwrap();
    //     if self.stars[&v].written == false {
    //         io::stdout().write_all(&format!("v {} {} {}\n", vpt[0], vpt[1], vpt[2]).as_bytes())?;
    //         self.stars.get_mut(&v).unwrap().smaid = self.smacount;
    //         self.smacount += 1;
    //         self.stars.get_mut(&v).unwrap().written = true;
    //     }
    //     let mut adjs: Vec<usize> = Vec::new();
    //     for each in self.stars[&v].link.iter() {
    //         adjs.push(*each);
    //     }
    //     for each in adjs {
    //         if (self.vertex_exists(each)) && (self.stars[&each].written == false) {
    //             let eachpt = self.get_point(each).unwrap();
    //             io::stdout().write_all(
    //                 &format!("v {} {} {}\n", eachpt[0], eachpt[1], eachpt[2]).as_bytes(),
    //             )?;
    //             self.stars.get_mut(&each).unwrap().smaid = self.smacount;
    //             self.smacount += 1;
    //             self.stars.get_mut(&each).unwrap().written = true;
    //         }
    //     }
    //     //-- write the faces/triangles
    //     for (i, each) in self.stars[&v].link.iter().enumerate() {
    //         let j = self.stars[&v].link.next_index(i);
    //         let smaj = self.stars[&v].link[j];
    //         if (*each != 0) //-- do not write out infinite triangles
    //             && (smaj != 0)
    //             && (self.vertex_exists(*each))
    //             && (self.vertex_exists(smaj))
    //         {
    //             io::stdout().write_all(
    //                 &format!(
    //                     "f {} {} {}\n",
    //                     self.stars[&v].smaid, self.stars[each].smaid, self.stars[&smaj].smaid
    //                 )
    //                 .as_bytes(),
    //             )?;
    //         }
    //     }
    //     io::stdout().write_all(&format!("x {}\n", self.stars[&v].smaid).as_bytes())?;
    //     Ok(())
    // }

    fn finalise_qt_root(&mut self) -> io::Result<()> {
        info!("Finalise the quadtree root cell");
        let q2: Vec<u8> = Vec::new();
        for ti in &self.qt.cells[&q2].ts {
            if self.is_triangle_finite(*ti) == false {
                continue;
            }
            io::stdout().write_all(
                &format!(
                    "f {} {} {}\n",
                    self.sma_ids[&self.ts[*ti][0]],
                    self.sma_ids[&self.ts[*ti][1]],
                    self.sma_ids[&self.ts[*ti][2]],
                )
                .as_bytes(),
            )?;
        }
        for vi in &self.qt.cells[&q2].pts {
            io::stdout().write_all(&format!("x {}\n", self.sma_ids[vi]).as_bytes())?;
        }
        Ok(())
    }

    // pub fn finalise_leftover_triangles(&mut self) -> io::Result<()> {
    //     let mut total: usize = 0;
    //     for i in &self.qt.gpts {
    //         for j in i {
    //             total += j.len();
    //         }
    //     }
    //     info!("Writing the {} vertices left in the DT", total);
    //     info!("DT # points: {}", self.number_of_vertices());

    //     //-- write the leftovers
    //     for i in 0..self.qt.griddim {
    //         for j in 0..self.qt.griddim {
    //             let vs = self.qt.gpts[i][j].clone();
    //             for v in vs {
    //                 match self.outputmode {
    //                     // write triangles
    //                     Outputmode::Sma => {
    //                         self.write_sma_one_vertex(v)
    //                             .expect("Failed to write vertex");
    //                         //-- flush it from QT and the DS
    //                         self.flush_star(v);
    //                     }

    //                     Outputmode::Stars => {
    //                         let re = self.adjacent_vertices_to_vertex(v).unwrap();
    //                         io::stdout().write_all(&format!("x {} {:?}\n", v, re).as_bytes())?;
    //                         //-- flush it from QT and the DS
    //                         self.flush_star(v);
    //                     }

    //                     //-- finalise the vertex, output both vertex info and star data
    //                     Outputmode::Both => {
    //                         self.write_stars_one_vertex(v)
    //                             .expect("Failed to write vertex");
    //                         //-- flush it from QT and the DS
    //                         self.flush_star(v);
    //                     }

    //                     _ => {
    //                         println!("_ output");
    //                     }
    //                 }
    //             }
    //         }
    //     }

    //     Ok(())
    // }

    // fn flush_star(&mut self, v: usize) -> bool {
    //     // self.stars.get_mut(&v).unwrap().active = false;
    //     // true
    //     let re = self.stars.remove(&v);
    //     if re.is_some() {
    //         true
    //     } else {
    //         println!("=== OH NO ===");
    //         false
    //     }
    // }

    // pub fn set_outputmode(&mut self, outmode: Outputmode) {
    //     self.outputmode = outmode;
    // }

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
        //-- checking if within snap_tol of one vertex, then no insert
        for i in 0..3 {
            if geom::distance2d_squared(&self.vs[self.ts[ti][i]], &p)
                <= (self.snaptol * self.snaptol)
            {
                return Err(self.ts[ti][i]);
            }
        }
        //-- ok we now insert the point in the data structure
        let pi: usize = self.vs.len();
        self.vs.push([px, py, pz]);

        //-- flip13()
        // self.print_ds();
        let incidents = self.flip13(pi, ti);
        // self.print_ds();

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
        let mut newcurt = self.ts.len() - 1;
        if self.is_triangle_finite(newcurt) == false {
            let it = self.ts[newcurt];
            let k = &it[0..3].iter().position(|&x| x == 0).unwrap();
            newcurt = it[k + 3];
        }
        self.curt = newcurt;
        Ok(pi)
    }

    fn flip13(&mut self, vi: usize, ti: usize) -> (usize, usize, usize) {
        let newi = self.ts.len();
        let t = self.ts[ti];
        let t0 = [vi, t[0], t[1], t[5], newi, newi + 1];
        let t1 = [vi, t[1], t[2], t[3], newi + 1, ti];
        let t2 = [vi, t[2], t[0], t[4], ti, newi];
        //-- update neighbours
        let mut n = t[3];
        if n != 0 {
            let a = self.ts[n];
            let a1 = &a[3..6].iter().position(|&x| x == ti).unwrap();
            if let Some(x) = self.ts.get_mut(n) {
                (*x)[*a1 + 3] = newi;
            }
        }
        n = t[4];
        if n != 0 {
            let a = self.ts[n];
            let a1 = &a[3..6].iter().position(|&x| x == ti).unwrap();
            if let Some(x) = self.ts.get_mut(n) {
                (*x)[*a1 + 3] = newi + 1;
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
        self.ts.insert(newi, t1);
        self.ts.insert(newi + 1, t2);
        (ti, newi, newi + 1)
    }

    // /// Returns the coordinates of the vertex v in a Vec [x,y,z]
    // pub fn get_point(&self, v: usize) -> Option<Vec<f64>> {
    //     if self.vertex_exists(v) == false {
    //         None
    //     } else {
    //         Some(self.stars[&v].pt.to_vec())
    //     }
    // }

    /// Returns a list (`Vec<usize>`) (ordered CCW) of the adjacent vertices.
    /// [`None`] if the vertex is not part of the triangulation.
    // pub fn adjacent_vertices_to_vertex(&self, v: usize) -> Option<Vec<usize>> {
    //     if self.vertex_exists(v) == false {
    //         return None;
    //     }
    //     let mut adjs: Vec<usize> = Vec::new();
    //     for each in self.stars[&v].link.iter() {
    //         adjs.push(*each);
    //     }
    //     Some(adjs)
    // }

    pub fn incident_triangles_to_vertex(&self, vi: usize) -> Option<Vec<usize>> {
        let mut l: Vec<usize> = Vec::new();
        let t0 = self.walk(&[self.vs[vi][0], self.vs[vi][1]]);
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
        self.vs.len() - 1
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

    // /// Returns the convex hull of the dataset, oriented CCW.
    // /// It is a list of vertex indices (first != last)
    // pub fn convex_hull(&self) -> Vec<usize> {
    //     let mut re: Vec<usize> = Vec::new();
    //     for x in self.stars[&0].link.iter() {
    //         re.push(*x);
    //     }
    //     re.reverse();
    //     re
    // }

    // /// Returns the size (ie the number of vertices) of the convex hull of the dataset
    // pub fn number_of_vertices_on_convex_hull(&self) -> usize {
    //     //-- number of finite vertices on the boundary of the convex hull
    //     if self.is_init == false {
    //         return 0;
    //     }
    //     return self.stars[&0].link.len();
    // }

    // /// Returns true if the vertex v is part of the boundary of the convex
    // /// hull of the dataset. False otherwise.
    // pub fn is_vertex_convex_hull(&self, v: usize) -> bool {
    //     if v == 0 {
    //         return false;
    //     }
    //     if self.vertex_exists(v) == false {
    //         return false;
    //     }
    //     self.stars[&v].link.contains_infinite_vertex()
    // }

    fn walk(&self, x: &[f64]) -> usize {
        //-- find a starting tr
        let tr = self.curt;

        //-- TODO: maybe based on QT size if larger than go to cell first?

        //-- 1. try walk from latest
        let re = self.walk_safe(x, tr);
        if re.is_some() {
            return re.unwrap();
        }

        //-- 2. try walk from one in the same cell
        // TODO: try walk from one in the same cell BROKEN
        // warn!("attempt to find one vertex in the grid cell and start from it");
        // let qtc = self.qt.get_qtc_from_xy(x[0], x[1]);
        // if self.qt.cells[&qtc].number_pts() > 0 {
        //     let v0 = self.qt.cells[&qtc].pts.iter().next().unwrap();
        //     // let v0 = self.qt.gpts[g.0][g.1].get(a[0]);
        //     let re = self.walk_safe(x, *v0);
        //     if re.is_some() {
        //         return re.unwrap();
        //     }
        // }

        //-- 3. try brute-force
        //-- TODO: this brute-force too?
        // let re2 = self.walk_bruteforce_closest_vertex_then_walksafe(x);
        // if re2.is_some() {
        //     return re2.unwrap();
        // }
        let re3 = self.walk_bruteforce_triangles(x);
        if re3.is_some() {
            return re3.unwrap();
        }

        //-- 4. we are outside the CH of the current dataset
        // warn!("point is outside the CH, finding closest point on the CH");
        let re4 = self.walk_bruteforce_outside_convex_hull(x);
        if re4.is_some() {
            return re4.unwrap();
        } else {
            error!("WALK FAILED MISERABLY :'(");
        }

        return 0;
    }

    fn walk_safe(&self, x: &[f64], mut tr: usize) -> Option<usize> {
        loop {
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
        for (id, t) in self.ts.iter().enumerate() {
            if self.is_triangle_finite(id) == true {
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
