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

/// A Triangle is a triplet of indices
pub struct Triangle {
    pub v: [usize; 3],
}

impl Triangle {
    /// Checks whether a Triangle is "infinite",
    /// ie if one its vertices is the infinite vertex
    fn is_infinite(&self) -> bool {
        if self.v[0] == 0 || self.v[1] == 0 || self.v[2] == 0 {
            return true;
        }
        return false;
    }
}

impl fmt::Display for Triangle {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({}, {}, {})", self.v[0], self.v[1], self.v[2])
    }
}

//----------------------
struct Link(Vec<usize>);

impl Link {
    fn new() -> Link {
        // Link(Vec::new())
        Link(Vec::with_capacity(8))
    }
    fn len(&self) -> usize {
        self.0.len()
    }
    fn is_empty(&self) -> bool {
        if self.0.len() == 0 {
            true
        } else {
            false
        }
    }
    fn add(&mut self, v: usize) {
        self.0.push(v);
    }
    fn insert_after_v(&mut self, v: usize, after: usize) {
        let pos = self.0.iter().position(|&x| x == after).unwrap();
        self.0.insert(pos + 1, v);
    }
    fn delete(&mut self, v: usize) {
        let re = self.0.iter().position(|&x| x == v);
        if re != None {
            self.0.remove(re.unwrap());
        }
    }
    fn replace(&mut self, v: usize, newv: usize) {
        let re = self.0.iter().position(|&x| x == v);
        if re != None {
            self.0[re.unwrap()] = newv;
            // self.0.remove(re.unwrap());
        }
    }
    fn infinite_first(&mut self) {
        let re = self.0.iter().position(|&x| x == 0);
        if re != None {
            let posinf = re.unwrap();
            if posinf == 0 {
                return;
            }
            let mut newstar: Vec<usize> = Vec::new();
            for j in posinf..self.0.len() {
                newstar.push(self.0[j]);
            }
            for j in 0..posinf {
                newstar.push(self.0[j]);
            }
            // println!("newstar: {:?}", newstar);
            self.0 = newstar;
        }
    }
    fn clear(&mut self) {
        self.0.clear();
    }

    fn contains_infinite_vertex(&self) -> bool {
        let pos = self.0.iter().position(|&x| x == 0);
        if pos == None {
            return false;
        } else {
            return true;
        }
    }

    fn next_index(&self, i: usize) -> usize {
        if i == (self.0.len() - 1) {
            0
        } else {
            i + 1
        }
    }

    fn prev_index(&self, i: usize) -> usize {
        if i == 0 {
            (self.0.len() - 1)
        } else {
            (i - 1)
        }
    }

    fn get_index(&self, v: usize) -> Option<usize> {
        return self.0.iter().position(|&x| x == v);
    }

    fn get_next_vertex(&self, v: usize) -> Option<usize> {
        let re = self.get_index(v);
        if re.is_none() {
            return None;
        }
        let pos = re.unwrap();
        if pos == (self.0.len() - 1) {
            return Some(self.0[0]);
        } else {
            return Some(self.0[(pos + 1)]);
        }
    }

    fn get_prev_vertex(&self, v: usize) -> Option<usize> {
        let re = self.get_index(v);
        if re.is_none() {
            return None;
        }
        let pos = re.unwrap();
        if pos == 0 {
            return Some(self.0[(self.0.len() - 1)]);
        } else {
            return Some(self.0[(pos - 1)]);
        }
    }

    fn iter(&self) -> Iter {
        Iter(Box::new(self.0.iter()))
    }
}

//-- taken from https://stackoverflow.com/questions/40668074/am-i-incorrectly-implementing-intoiterator-for-a-reference-or-is-this-a-rust-bug
struct Iter<'a>(Box<dyn Iterator<Item = &'a usize> + 'a>);

impl<'a> Iterator for Iter<'a> {
    type Item = &'a usize;

    fn next(&mut self) -> Option<Self::Item> {
        self.0.next()
    }
}

impl std::ops::Index<usize> for Link {
    type Output = usize;
    fn index(&self, idx: usize) -> &usize {
        &self.0[idx as usize]
    }
}

impl fmt::Display for Link {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        // fmt.write_str("pt: {}\n", self.pt)?;
        fmt.write_str(&format!("link: {:?}\n", self.0))?;
        Ok(())
    }
}

/// A triangulation is a collection of Stars, each Star has its (x,y,z)
/// and a Link (an array of adjacent vertices, ordered CCW)
struct Star {
    pub pt: [f64; 3],
    pub link: Link,
    pub active: bool,
}

impl Star {
    pub fn new(x: f64, y: f64, z: f64) -> Star {
        let l = Link::new();
        Star {
            pt: [x, y, z],
            link: l,
            active: true,
        }
    }
}

//----------------------
struct Quadtree {
    pub gpts: Vec<Vec<HashSet<usize>>>,
    pub gfinal: HashSet<Vec<u8>>,
    pub minx: f64,
    pub miny: f64,
    pub cellsize: usize,
    pub griddim: usize,
    depth: u32,
}

impl Quadtree {
    pub fn new() -> Quadtree {
        let g: Vec<Vec<HashSet<usize>>> = Vec::new();
        let gf: HashSet<Vec<u8>> = HashSet::new();
        Quadtree {
            gpts: g,
            gfinal: gf,
            cellsize: 0,
            griddim: 0,
            minx: std::f64::MAX,
            miny: std::f64::MAX,
            depth: 0,
        }
    }

    pub fn init(&mut self, griddim: usize) {
        self.griddim = griddim;
        self.depth = (griddim as f64).log(2.0).ceil() as u32;
        self.gpts.resize(griddim, Vec::new());
        for each in &mut self.gpts {
            each.resize(griddim, HashSet::new());
        }
    }

    pub fn get_cell_count(&self, gx: usize, gy: usize) -> Option<usize> {
        if gx >= self.griddim || gy > self.griddim {
            None
        } else {
            Some(self.gpts[gx][gy].len())
        }
    }

    pub fn get_cell_pts(&self, gx: usize, gy: usize) -> Vec<usize> {
        let mut l: Vec<usize> = Vec::new();
        for each in self.gpts[gx][gy].iter() {
            l.push(*each);
        }
        l
    }

    pub fn get_cell_pts_qtc(&self, qtc: &Vec<u8>) -> Vec<usize> {
        let (gx, gy) = self.qtc2gxgy(qtc);
        let mut l: Vec<usize> = Vec::new();
        for each in self.gpts[gx][gy].iter() {
            l.push(*each);
        }
        l
    }

    pub fn get_all_gxgy_from_qtc(&self, qtc: &Vec<u8>) -> Vec<(usize, usize)> {
        let mut re: Vec<(usize, usize)> = Vec::new();
        if qtc.len() as u32 == self.depth {
            re.push(self.qtc2gxgy(qtc));
            return re;
        }
        let mut stack: Vec<Vec<u8>> = Vec::new();
        stack.push(qtc.to_vec());
        while stack.is_empty() == false {
            let q = stack.pop().unwrap();
            if q.len() as u32 == self.depth {
                re.push(self.qtc2gxgy(&q));
            } else {
                for i in 0..4 {
                    let mut q2 = q.to_vec();
                    q2.push(i);
                    stack.push(q2);
                }
            }
        }
        re
    }

    /// Returns the qtc of the cell/s that is/are finalised, it can be only the
    /// current one or its parent (or parent of that one; check recursively)
    pub fn finalise_cell(&mut self, gx: usize, gy: usize) -> Vec<u8> {
        let qtc = self.get_cell_qtc(gx, gy);
        let mut q2 = vec![0; qtc.len()];
        q2.clone_from_slice(&qtc);
        //-- add to the hashset
        self.gfinal.insert(qtc);
        //-- check parents recursively and add them to gfinals
        let mut done = false;
        while !done {
            if q2.is_empty() == true {
                // done = true;
                break;
            }
            if self.finalise_parent(&q2) == true {
                q2.pop();
            } else {
                done = true;
            }
        }
        println!("{:?}", q2);
        q2
    }

    fn finalise_parent(&mut self, qtc: &Vec<u8>) -> bool {
        let mut q2 = vec![0; qtc.len() - 1];
        q2.clone_from_slice(&qtc[..qtc.len() - 1]);
        if self.gfinal.contains(&q2) == true {
            return true;
        }
        let f = true;
        q2.push(0);
        if self.gfinal.contains(&q2) == false {
            return false;
        }
        q2.pop();
        q2.push(1);
        if self.gfinal.contains(&q2) == false {
            return false;
        }
        q2.pop();
        q2.push(2);
        if self.gfinal.contains(&q2) == false {
            return false;
        }
        q2.pop();
        q2.push(3);
        if self.gfinal.contains(&q2) == false {
            return false;
        }
        q2.pop();
        self.gfinal.insert(q2);
        true
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

    fn get_cell_gxgy(&self, x: f64, y: f64) -> (usize, usize) {
        (
            ((x - self.minx) / self.cellsize as f64) as usize,
            ((y - self.miny) / self.cellsize as f64) as usize,
        )
    }

    fn get_cell_qtc(&self, gx: usize, gy: usize) -> Vec<u8> {
        // println!("{:0>4b}", gx);
        // println!("{}", gy.leading_zeros());
        // println!("{:0>4b}", gx & 0b0100);
        // println!("{}", gx & 0b0100);
        let mut mask: usize = 2_usize.pow(self.depth - 1);
        // println!("mask: {}", mask);
        // println!("mask: {:0>4b}", mask);
        // let mut a: usize = 0b1000;
        let mut re: Vec<u8> = Vec::new();
        for i in 0..self.depth {
            // println!("i:{} -- {:0>4b}", i, mask);
            // println!("{}", gx & mask == 2_usize.pow(self.depth - i - 1));
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

    fn is_cell_final(&self, gx: usize, gy: usize) -> bool {
        self.gfinal.contains(&self.get_cell_qtc(gx, gy))
        // TODO: delete the leaves that are delete by their parents?
    }

    fn is_cell_final_qtc(&self, qtc: &Vec<u8>) -> bool {
        self.gfinal.contains(qtc)
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
        if c[0] == 0 {
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
    stars: HashMap<usize, Star>,
    qt: Quadtree,
    snaptol: f64,
    cur: usize,
    is_init: bool,
    robust_predicates: bool,
    theid: usize, //-- to assign id to new vertices, since some are flushed
}

impl Triangulation {
    //-- new
    pub fn new() -> Triangulation {
        let mut s: HashMap<usize, Star> = HashMap::new();
        s.insert(0, Star::new(-99999.99999, -99999.99999, -99999.99999));
        let q = Quadtree::new();
        Triangulation {
            qt: q,
            theid: 1,
            stars: s,
            snaptol: 0.001,
            cur: 0,
            is_init: false,
            robust_predicates: true,
        }
    }

    pub fn finalise_cell(&mut self, gx: usize, gy: usize) -> io::Result<()> {
        info!(
            "Cell {}--{} finalised ({} vertices)",
            gx,
            gy,
            self.qt.get_cell_count(gx, gy).unwrap()
        );

        let qtc = self.qt.finalise_cell(gx, gy);

        if qtc.is_empty() == true {
            return Ok(());
        }
        //-- nothing to do if single cell finalised and it's empty
        if (qtc.len() as u32 == self.qt.depth) && (self.qt.get_cell_count(gx, gy).unwrap() == 0) {
            return Ok(());
        }

        let allfcells = self.qt.get_all_gxgy_from_qtc(&qtc);

        let mut gbbox: [f64; 4] = [0.0, 0.0, 0.0, 0.0];
        self.qt.get_cell_bbox_qtc(&qtc, &mut gbbox);
        for c in allfcells {
            let mut finpts: HashSet<usize> = HashSet::new();
            for theid in self.qt.gpts[c.0][c.1].iter() {
                let re = self.adjacent_vertices_to_vertex(*theid).unwrap();
                let mut fin: bool = true;
                for v in re {
                    // if self.qt.gpts[gx][gy].contains(&v) == false {
                    let p = self.get_point(v);
                    if p.is_none() {
                        continue;
                    }
                    // if p.is_none()
                    if geom::point_in_box(&self.get_point(v).unwrap(), &gbbox) == false {
                        fin = false;
                        break;
                    }
                }
                if fin == true {
                    // println!("=>{}", theid);
                    //-- check every triangle for encroachment
                    //-- TODO: should incident function only returns non-final triangle?
                    let lts = self.incident_triangles_to_vertex(*theid).unwrap();
                    for t in &lts {
                        // println!("t {}", t);
                        if self.is_vertex_convex_hull(t.v[1]) && self.is_vertex_convex_hull(t.v[2])
                        {
                            fin = false;
                            break;
                        }
                        let a = self.get_point(t.v[0]);
                        let b = self.get_point(t.v[1]);
                        let c = self.get_point(t.v[2]);
                        if a.is_none() == false
                            && b.is_none() == false
                            && c.is_none() == false
                            && geom::circumcentre_encroach_bbox(
                                &a.unwrap(),
                                &b.unwrap(),
                                &c.unwrap(),
                                &gbbox,
                            ) == true
                        {
                            fin = false;
                            break;
                        }
                    }
                }
                if fin == true {
                    finpts.insert(*theid);
                }
            }
            for each in &finpts {
                let p = self.get_point(*each).unwrap();
                io::stdout().write_all(
                    &format!(
                        "v {} {} {} {} {:?}\n",
                        *each,
                        p[0],
                        p[1],
                        p[2],
                        self.adjacent_vertices_to_vertex(*each).unwrap()
                    )
                    .as_bytes(),
                )?;
                self.qt.gpts[gx][gy].remove(each);
            }
            for each in &finpts {
                self.flush_star(*each);
            }
        }
        Ok(())
    }

    pub fn finalise_leftover_triangles(&mut self) -> io::Result<()> {
        let mut total: usize = 0;
        for w in &self.qt.gpts {
            for h in w {
                total += h.len();
            }
        }
        info!("Writing the {} vertices left in the DT", total);
        info!("DT # points: {}", self.number_of_vertices());
        //-- write the leftovers
        for w in &self.qt.gpts {
            for h in w {
                for each in h.iter() {
                    let p = self.get_point(*each).unwrap();
                    io::stdout().write_all(
                        &format!(
                            "v {} {} {} {} {:?}\n",
                            *each,
                            p[0],
                            p[1],
                            p[2],
                            self.adjacent_vertices_to_vertex(*each).unwrap()
                        )
                        .as_bytes(),
                    )?;
                }
            }
        }
        Ok(())
    }

    fn flush_star(&mut self, v: usize) -> bool {
        // self.stars.get_mut(&v).unwrap().active = false;
        // true
        let re = self.stars.remove(&v);
        if re.is_some() {
            true
        } else {
            println!("=== OH NO ===");
            false
        }
    }

    pub fn set_cellsize(&mut self, c: usize) {
        self.qt.cellsize = c;
    }

    pub fn set_bbox(&mut self, mx: f64, my: f64) {
        self.qt.minx = mx;
        self.qt.miny = my;
    }

    pub fn set_grid_dimensions(&mut self, s: usize) {
        self.qt.init(s);
    }

    fn insert_one_pt_init_phase(&mut self, x: f64, y: f64, z: f64) -> Result<usize, usize> {
        let p: [f64; 3] = [x, y, z];
        for i in 1..self.stars.len() {
            if geom::distance2d_squared(&self.stars[&i].pt, &p) <= (self.snaptol * self.snaptol) {
                return Err(i);
            }
        }
        //-- add point to Triangulation and create its empty star
        self.stars.insert(self.theid, Star::new(x, y, z));
        self.theid += 1;
        //-- form the first triangles (finite + infinite)
        let l = self.stars.len();
        if l >= 4 {
            let a = l - 3;
            let b = l - 2;
            let c = l - 1;
            let re = geom::orient2d(
                &self.stars[&a].pt,
                &self.stars[&b].pt,
                &self.stars[&c].pt,
                self.robust_predicates,
            );
            if re == 1 {
                // println!("init: ({},{},{})", a, b, c);
                self.stars.get_mut(&0).unwrap().link.add(a);
                self.stars.get_mut(&0).unwrap().link.add(c);
                self.stars.get_mut(&0).unwrap().link.add(b);
                self.stars.get_mut(&a).unwrap().link.add(0);
                self.stars.get_mut(&a).unwrap().link.add(b);
                self.stars.get_mut(&a).unwrap().link.add(c);
                self.stars.get_mut(&b).unwrap().link.add(0);
                self.stars.get_mut(&b).unwrap().link.add(c);
                self.stars.get_mut(&b).unwrap().link.add(a);
                self.stars.get_mut(&c).unwrap().link.add(0);
                self.stars.get_mut(&c).unwrap().link.add(a);
                self.stars.get_mut(&c).unwrap().link.add(b);
                self.is_init = true;
            } else if re == -1 {
                // println!("init: ({},{},{})", a, c, b);
                self.stars.get_mut(&0).unwrap().link.add(a);
                self.stars.get_mut(&0).unwrap().link.add(b);
                self.stars.get_mut(&0).unwrap().link.add(c);
                self.stars.get_mut(&a).unwrap().link.add(0);
                self.stars.get_mut(&a).unwrap().link.add(c);
                self.stars.get_mut(&a).unwrap().link.add(b);
                self.stars.get_mut(&b).unwrap().link.add(0);
                self.stars.get_mut(&b).unwrap().link.add(a);
                self.stars.get_mut(&b).unwrap().link.add(c);
                self.stars.get_mut(&c).unwrap().link.add(0);
                self.stars.get_mut(&c).unwrap().link.add(b);
                self.stars.get_mut(&c).unwrap().link.add(a);
                self.is_init = true;
            }
        }
        self.cur = l - 1;
        if self.is_init == true {
            //-- insert the previous vertices in the dt
            for j in 1..(l - 3) {
                let tr = self.walk(&self.stars[&j].pt);
                // println!("found tr: {}", tr);
                self.flip13(j, &tr);
                self.update_dt(j);
            }
        }
        Ok(self.cur)
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

    pub fn insert(&mut self, pts: &Vec<Vec<f64>>) {
        let mut duplicates = 0;
        for each in pts {
            if (each.len() < 2) || (each.len() > 3) {
                panic!(
                    "Point {:?} should be 2D or 3D (and is now {}D).",
                    each,
                    each.len()
                );
            } else {
                let re;
                if each.len() == 2 {
                    re = self.insert_one_pt(each[0], each[1], 0.0);
                } else {
                    re = self.insert_one_pt(each[0], each[1], each[2]);
                }
                match re {
                    Ok(_x) => continue,
                    Err(_e) => duplicates = duplicates + 1,
                }
            }
        }
    }

    pub fn insert_one_pt_with_grid(&mut self, px: f64, py: f64, pz: f64) -> Result<usize, usize> {
        let re = self.insert_one_pt(px, py, pz);
        if re.is_ok() {
            let x = re.unwrap();
            let g = self.qt.get_cell_gxgy(px, py);
            self.qt.gpts[g.0][g.1].insert(x);
        };
        re
    }

    //-- insert_one_pt
    fn insert_one_pt(&mut self, px: f64, py: f64, pz: f64) -> Result<usize, usize> {
        // println!("-->{}", p);
        if self.is_init == false {
            return self.insert_one_pt_init_phase(px, py, pz);
        }
        //-- walk
        let p: [f64; 3] = [px, py, pz];
        let tr = self.walk(&p);
        // println!("STARTING TR: {}", tr);
        if geom::distance2d_squared(&self.stars[&tr.v[0]].pt, &p) <= (self.snaptol * self.snaptol) {
            return Err(tr.v[0]);
        }
        if geom::distance2d_squared(&self.stars[&tr.v[1]].pt, &p) <= (self.snaptol * self.snaptol) {
            return Err(tr.v[1]);
        }
        if geom::distance2d_squared(&self.stars[&tr.v[2]].pt, &p) <= (self.snaptol * self.snaptol) {
            return Err(tr.v[2]);
        }
        //-- ok we now insert the point in the data structure
        //-- TODO: remove this for delete when hash is used
        let pi: usize = self.theid;
        self.stars.insert(self.theid, Star::new(px, py, pz));
        self.theid += 1;
        //-- flip13()
        self.flip13(pi, &tr);
        //-- update_dt()
        self.update_dt(pi);
        self.cur = pi;
        Ok(pi)
    }

    fn update_dt(&mut self, pi: usize) {
        // println!("--> Update DT");
        let mut mystack: Vec<Triangle> = Vec::new();
        let l = &self.stars.get_mut(&pi).unwrap().link;
        mystack.push(Triangle {
            v: [pi, l[0], l[1]],
        });
        mystack.push(Triangle {
            v: [pi, l[1], l[2]],
        });
        mystack.push(Triangle {
            v: [pi, l[2], l[0]],
        });

        loop {
            let tr = match mystack.pop() {
                None => break,
                Some(x) => x,
            };
            let opposite = self.get_opposite_vertex(&tr);
            //-- TODO: danger danger
            if self.stars.contains_key(&opposite) == false {
                continue;
            }
            //----------------------
            // println!("stacked: {} {}", tr, opposite);

            if tr.is_infinite() == true {
                let mut a: i8 = 0;
                if tr.v[0] == 0 {
                    a = geom::orient2d(
                        &self.stars[&opposite].pt,
                        &self.stars[&tr.v[1]].pt,
                        &self.stars[&tr.v[2]].pt,
                        self.robust_predicates,
                    );
                } else if tr.v[1] == 0 {
                    a = geom::orient2d(
                        &self.stars[&tr.v[0]].pt,
                        &self.stars[&opposite].pt,
                        &self.stars[&tr.v[2]].pt,
                        self.robust_predicates,
                    );
                } else if tr.v[2] == 0 {
                    a = geom::orient2d(
                        &self.stars[&tr.v[0]].pt,
                        &self.stars[&tr.v[1]].pt,
                        &self.stars[&opposite].pt,
                        self.robust_predicates,
                    );
                }
                // println!("TODO: INCIRCLE FOR INFINITY {}", a);
                if a > 0 {
                    // println!("FLIPPED0 {} {}", tr, opposite);
                    let (ret0, ret1) = self.flip22(&tr, opposite);
                    mystack.push(ret0);
                    mystack.push(ret1);
                }
            } else {
                if opposite == 0 {
                    //- if insertion on CH then break the edge, otherwise do nothing
                    //-- TODO sure the flips are okay here?
                    if geom::orient2d(
                        &self.stars[&tr.v[0]].pt,
                        &self.stars[&tr.v[1]].pt,
                        &self.stars[&tr.v[2]].pt,
                        self.robust_predicates,
                    ) == 0
                    {
                        // println!("FLIPPED1 {} {}", tr, 0);
                        let (ret0, ret1) = self.flip22(&tr, 0);
                        mystack.push(ret0);
                        mystack.push(ret1);
                    }
                } else {
                    if geom::incircle(
                        &self.stars[&tr.v[0]].pt,
                        &self.stars[&tr.v[1]].pt,
                        &self.stars[&tr.v[2]].pt,
                        &self.stars[&opposite].pt,
                        self.robust_predicates,
                    ) > 0
                    {
                        // println!("FLIPPED2 {} {}", tr, opposite);
                        let (ret0, ret1) = self.flip22(&tr, opposite);
                        mystack.push(ret0);
                        mystack.push(ret1);
                    }
                }
            }
        }
    }

    fn flip13(&mut self, pi: usize, tr: &Triangle) {
        let l = &mut self.stars.get_mut(&pi).unwrap().link;
        l.add(tr.v[0]);
        l.add(tr.v[1]);
        l.add(tr.v[2]);
        self.stars
            .get_mut(&tr.v[0])
            .unwrap()
            .link
            .insert_after_v(pi, tr.v[1]);
        self.stars
            .get_mut(&tr.v[1])
            .unwrap()
            .link
            .insert_after_v(pi, tr.v[2]);
        self.stars
            .get_mut(&tr.v[2])
            .unwrap()
            .link
            .insert_after_v(pi, tr.v[0]);
        //-- put infinite vertex first in list
        // self.stars[pi].link.infinite_first();
    }

    /// Returns the coordinates of the vertex v in a Vec [x,y,z]
    pub fn get_point(&self, v: usize) -> Option<Vec<f64>> {
        if self.vertex_exists(v) == false {
            None
        } else {
            Some(self.stars[&v].pt.to_vec())
        }
    }

    pub fn adjacent_triangles_to_triangle(&self, tr: &Triangle) -> Option<Vec<Triangle>> {
        if self.is_triangle(&tr) == false || tr.is_infinite() == true {
            return None;
        }
        let mut trs: Vec<Triangle> = Vec::new();
        let mut opp = self.stars[&tr.v[2]].link.get_next_vertex(tr.v[1]).unwrap();
        if opp != 0 {
            trs.push(Triangle {
                v: [tr.v[1], opp, tr.v[2]],
            });
        }
        opp = self.stars[&tr.v[0]].link.get_next_vertex(tr.v[2]).unwrap();
        if opp != 0 {
            trs.push(Triangle {
                v: [tr.v[2], opp, tr.v[0]],
            });
        }
        opp = self.stars[&tr.v[1]].link.get_next_vertex(tr.v[0]).unwrap();
        if opp != 0 {
            trs.push(Triangle {
                v: [tr.v[0], opp, tr.v[1]],
            });
        }
        Some(trs)
    }

    /// Returns a Vec of Triangles (finite + infinite) to the vertex v.
    /// If v doesn't exist, then [`None`] is returned.
    pub fn incident_triangles_to_vertex(&self, v: usize) -> Option<Vec<Triangle>> {
        if self.vertex_exists(v) == false {
            return None;
        }
        let mut trs: Vec<Triangle> = Vec::new();
        for (i, each) in self.stars[&v].link.iter().enumerate() {
            let j = self.stars[&v].link.next_index(i);
            trs.push(Triangle {
                v: [v, *each, self.stars[&v].link[j]],
            });
        }
        Some(trs)
    }

    /// Returns the degree of a vertex, [`None`] is it doesn't exist.
    pub fn degree(&self, v: usize) -> Option<usize> {
        if self.vertex_exists(v) == false {
            return None;
        }
        Some(self.stars[&v].link.len())
    }

    /// Returns a list (`Vec<usize>`) (ordered CCW) of the adjacent vertices.
    /// [`None`] if the vertex is not part of the triangulation.
    pub fn adjacent_vertices_to_vertex(&self, v: usize) -> Option<Vec<usize>> {
        if self.vertex_exists(v) == false {
            return None;
        }
        let mut adjs: Vec<usize> = Vec::new();
        for each in self.stars[&v].link.iter() {
            adjs.push(*each);
        }
        Some(adjs)
    }

    /// Returns whether a triplet of indices is a Triangle in the triangulation.
    pub fn is_triangle(&self, tr: &Triangle) -> bool {
        // TODO: what about infinite triangles?
        let re = self.stars[&tr.v[0]].link.get_next_vertex(tr.v[1]);
        if re.is_none() {
            return false;
        } else {
            if re.unwrap() == tr.v[2] {
                return true;
            } else {
                return false;
            }
        }
    }

    pub fn statistics_degree(&self) -> (f64, usize, usize) {
        let mut total: f64 = 0.0;
        let mut min: usize = usize::max_value();
        let mut max: usize = usize::min_value();
        for i in 1..self.stars.len() {
            total = total + self.stars[&i].link.len() as f64;
            if self.stars[&i].link.len() > max {
                max = self.stars[&i].link.len();
            }
            if self.stars[&i].link.len() < min {
                min = self.stars[&i].link.len();
            }
        }
        total = total / (self.stars.len() - 2) as f64;
        return (total, min, max);
    }

    /// Returns number of finite vertices in the triangulation.
    pub fn number_of_vertices(&self) -> usize {
        //-- number of finite vertices
        (self.stars.len() - 1)
    }

    /// Returns number of finite triangles in the triangulation.
    pub fn number_of_triangles(&self) -> usize {
        //-- number of finite triangles
        let mut count: usize = 0;
        for (i, star) in &self.stars {
            for (j, value) in star.link.iter().enumerate() {
                if i < value {
                    let k = star.link[star.link.next_index(j)];
                    if i < &k {
                        let tr = Triangle { v: [*i, *value, k] };
                        if tr.is_infinite() == false {
                            count = count + 1;
                        }
                    }
                }
            }
        }
        count
    }

    /// Returns the convex hull of the dataset, oriented CCW.
    /// It is a list of vertex indices (first != last)
    pub fn convex_hull(&self) -> Vec<usize> {
        let mut re: Vec<usize> = Vec::new();
        for x in self.stars[&0].link.iter() {
            re.push(*x);
        }
        re.reverse();
        re
    }

    /// Returns the size (ie the number of vertices) of the convex hull of the dataset
    pub fn number_of_vertices_on_convex_hull(&self) -> usize {
        //-- number of finite vertices on the boundary of the convex hull
        if self.is_init == false {
            return 0;
        }
        return self.stars[&0].link.len();
    }

    /// Returns true if the vertex v is part of the boundary of the convex
    /// hull of the dataset. False otherwise.
    pub fn is_vertex_convex_hull(&self, v: usize) -> bool {
        if v == 0 {
            return false;
        }
        if self.vertex_exists(v) == false {
            return false;
        }
        self.stars[&v].link.contains_infinite_vertex()
    }

    /// Returns, if it exists, the Triangle containing (px,py).
    /// If it is direction on a vertex/edge, then one is randomly chosen.
    pub fn locate(&self, px: f64, py: f64) -> Option<Triangle> {
        let p: [f64; 3] = [px, py, 0.0];
        let re = self.walk(&p);
        match re.is_infinite() {
            true => None,
            false => Some(re),
        }
    }

    // Returns closest point (in 2D) to a query point (x,y).
    // if (x,y) is outside the convex hull [`None`]
    pub fn closest_point(&self, px: f64, py: f64) -> Option<usize> {
        let re = self.locate(px, py);
        if re.is_none() == true {
            return None;
        }
        let p: [f64; 3] = [px, py, 0.0];
        let tr = re.unwrap();
        let mut d = std::f64::MAX;
        let mut closest: usize = 0;
        //-- 1. find triangle and closest vertex from the 3
        for each in tr.v.iter() {
            // println!("{}", each);
            let dtmp = geom::distance2d_squared(&self.stars[each].pt, &p);
            if dtmp < d {
                d = dtmp;
                closest = *each;
            }
        }
        for each in self.stars[&closest].link.iter() {
            let dtmp = geom::distance2d_squared(&self.stars[each].pt, &p);
            if dtmp < d {
                d = dtmp;
                closest = *each;
            }
        }
        Some(closest)
    }

    fn walk(&self, x: &[f64]) -> Triangle {
        //-- find a starting tr
        let mut cur = self.cur;

        //-- 1. try walk from latest
        let re = self.walk_safe(x, cur);
        if re.is_some() {
            return re.unwrap();
        }
        warn!("'normal' walk() failed");

        //-- 2. try walk from one in the same cell
        warn!("attempt to find one vertex in the grid cell and start from it");
        let g = self.qt.get_cell_gxgy(x[0], x[1]);
        if self.qt.gpts[g.0][g.1].len() > 0 {
            cur = *self.qt.gpts[g.0][g.1].iter().next().unwrap();
            let re = self.walk_safe(x, cur);
            if re.is_some() {
                return re.unwrap();
            }
        }

        //-- 3. try brute-force on triangles
        warn!("attempt to find one triangle brute-force");
        let re2 = self.walk_bruteforce_triangles(x);
        if re2.is_some() {
            return re2.unwrap();
        }

        //-- 4. we are outside the CH of the current dataset
        warn!("point is outside the CH, finding closest point on the CH");
        let re4 = self.walk_bruteforce_vertices(x);
        // for key in self.stars.keys() {
        //     println!("{:?}", key);
        // }
        if re4.is_some() {
            return re4.unwrap();
        } else {
            error!("WALK FAILED MISERABLY :'(");
        }

        let tr = Triangle { v: [0, 0, 0] };
        return tr;
    }

    fn walk_safe(&self, x: &[f64], cur: usize) -> Option<Triangle> {
        let mut tr = Triangle { v: [0, 0, 0] };
        let re = &self.stars.get(&cur);
        if re.is_none() {
            return None;
        }

        //--
        // println!(
        //     "walk_safe: ({}) -- {}",
        //     cur,
        //     self.stars.get(&cur).unwrap().link
        // );
        // let a: Vec<&usize> = self.stars.keys().collect();
        // println!("{:?}", a);
        //--

        tr.v[0] = cur;
        let l = &re.unwrap().link;
        let mut b = false;
        for i in 0..(l.len() - 1) {
            if (l[i] != 0)
                && (l[i + 1] != 0)
                && (self.stars.contains_key(&l[i]) == true)
                && (self.stars.contains_key(&l[i + 1]) == true)
            {
                tr.v[1] = l[i];
                tr.v[2] = l[i + 1];
                b = true;
                break;
            }
        }
        if b == false {
            info!("Cannot find a starting finite triangle.");
            let pr: usize = l.get_prev_vertex(0).unwrap();
            tr.v[0] = cur;
            tr.v[1] = pr;
            tr.v[2] = 0;
            return Some(tr);
        }

        //-- 2. order it such that tr0-tr1-x is CCW
        if geom::orient2d(
            &self.stars[&tr.v[0]].pt,
            &self.stars[&tr.v[1]].pt,
            &x,
            self.robust_predicates,
        ) == -1
        {
            if geom::orient2d(
                &self.stars[&tr.v[1]].pt,
                &self.stars[&tr.v[2]].pt,
                &x,
                self.robust_predicates,
            ) != -1
            {
                let tmp: usize = tr.v[0];
                tr.v[0] = tr.v[1];
                tr.v[1] = tr.v[2];
                tr.v[2] = tmp;
            } else {
                let tmp: usize = tr.v[1];
                tr.v[0] = tr.v[2];
                tr.v[1] = tr.v[0];
                tr.v[2] = tmp;
            }
        }

        //-- 3. start the walk
        //-- we know that tr0-tr1-x is CCW
        loop {
            if tr.is_infinite() == true {
                break;
            }
            let a = &self.stars.get(&tr.v[0]);
            if a.is_none() {
                return None;
            }
            let b = &self.stars.get(&tr.v[1]);
            if b.is_none() {
                return None;
            }
            let c = &self.stars.get(&tr.v[2]);
            if c.is_none() {
                return None;
            }
            if geom::orient2d(
                &self.stars.get(&tr.v[1]).unwrap().pt,
                &self.stars.get(&tr.v[2]).unwrap().pt,
                &x,
                self.robust_predicates,
            ) != -1
            {
                if geom::orient2d(
                    &self.stars.get(&tr.v[2]).unwrap().pt,
                    &self.stars.get(&tr.v[0]).unwrap().pt,
                    &x,
                    self.robust_predicates,
                ) != -1
                {
                    break;
                } else {
                    //-- walk to incident to tr1,tr2
                    // println!("here");
                    let prev = self
                        .stars
                        .get(&tr.v[2])
                        .unwrap()
                        .link
                        .get_prev_vertex(tr.v[0])
                        .unwrap();
                    tr.v[1] = tr.v[2];
                    tr.v[2] = prev;
                }
            } else {
                //-- walk to incident to tr1,tr2
                // a.iter().position(|&x| x == 2), Some(1)
                let prev = self
                    .stars
                    .get(&tr.v[1])
                    .unwrap()
                    .link
                    .get_prev_vertex(tr.v[2])
                    .unwrap();
                tr.v[0] = tr.v[2];
                tr.v[2] = prev;
            }
        }
        return Some(tr);
    }

    fn walk_bruteforce_vertices(&self, x: &[f64]) -> Option<Triangle> {
        //-- find closest vertex that is on the CH
        let mut dmin: f64 = std::f64::MAX;
        let mut vmin: usize = 0;
        for i in self.stars.keys() {
            if (*i != 0) && (self.is_vertex_convex_hull(*i)) {
                let d = geom::distance2d_squared(x, &self.stars.get(i).unwrap().pt);
                if d < dmin {
                    dmin = d;
                    vmin = *i;
                }
            }
        }
        info!("brute-force ON CONVEX HULL");
        let mut tr = Triangle { v: [0, 0, 0] };
        let l = &self.stars.get(&vmin).unwrap().link;
        let pr: usize = l.get_prev_vertex(0).unwrap();
        tr.v[0] = vmin;
        tr.v[1] = pr;
        tr.v[2] = 0;
        return Some(tr);
        // self.walk_safe(x, vmin)
    }

    fn walk_bruteforce_triangles(&self, x: &[f64]) -> Option<Triangle> {
        for (i, star) in &self.stars {
            for (j, value) in star.link.iter().enumerate() {
                if (i < value) && (self.stars.contains_key(&value) == true) {
                    let k = star.link[star.link.next_index(j)];
                    if (i < &k) && (self.stars.contains_key(&k) == true) {
                        let tr = Triangle { v: [*i, *value, k] };
                        if tr.is_infinite() == false {
                            if geom::intriangle(
                                &self.stars.get(&tr.v[0]).unwrap().pt,
                                &self.stars.get(&tr.v[1]).unwrap().pt,
                                &self.stars.get(&tr.v[2]).unwrap().pt,
                                &x,
                                self.robust_predicates,
                            ) == 1
                            {
                                return Some(tr);
                            }
                        }
                    }
                }
            }
        }
        return None;
    }

    fn walk_old(&self, x: &[f64]) -> Triangle {
        let g = self.qt.get_cell_gxgy(x[0], x[1]);
        //-- find the starting tr
        let mut cur = self.cur;
        if !self.stars.contains_key(&cur) {
            error!("walk(): cur={} doesn't exist", cur);
            //-- fetch one from the cell
            if self.qt.gpts[g.0][g.1].len() > 0 {
                cur = *self.qt.gpts[g.0][g.1].iter().next().unwrap();
            } else {
                // cur = *self.stars.keys().next().unwrap();
                info!("Cell has no vertex in it.");
                cur = *self.stars.keys().next().unwrap();
            }
        }
        let mut tr = Triangle { v: [0, 0, 0] };
        // println!("cur: {}", cur);

        //-- 1. find a finite triangle
        tr.v[0] = cur;
        let l = &self.stars.get(&cur).unwrap().link;
        let mut b = false;
        for i in 0..(l.len() - 1) {
            if (l[i] != 0)
                && (l[i + 1] != 0)
                && (self.stars.contains_key(&l[i]) == true)
                && (self.stars.contains_key(&l[i + 1]) == true)
            {
                tr.v[1] = l[i];
                tr.v[2] = l[i + 1];
                b = true;
                break;
            }
        }
        if b == false {
            info!("Cannot find a starting triangle.");
        }

        //-- 2. order it such that tr0-tr1-x is CCW
        if geom::orient2d(
            &self.stars[&tr.v[0]].pt,
            &self.stars[&tr.v[1]].pt,
            &x,
            self.robust_predicates,
        ) == -1
        {
            if geom::orient2d(
                &self.stars[&tr.v[1]].pt,
                &self.stars[&tr.v[2]].pt,
                &x,
                self.robust_predicates,
            ) != -1
            {
                let tmp: usize = tr.v[0];
                tr.v[0] = tr.v[1];
                tr.v[1] = tr.v[2];
                tr.v[2] = tmp;
            } else {
                let tmp: usize = tr.v[1];
                tr.v[0] = tr.v[2];
                tr.v[1] = tr.v[0];
                tr.v[2] = tmp;
            }
        }

        info!("here");
        //-- 3. start the walk
        //-- we know that tr0-tr1-x is CCW
        loop {
            if tr.is_infinite() == true {
                break;
            }
            if geom::orient2d(
                &self.stars[&tr.v[1]].pt,
                &self.stars[&tr.v[2]].pt,
                &x,
                self.robust_predicates,
            ) != -1
            {
                if geom::orient2d(
                    &self.stars[&tr.v[2]].pt,
                    &self.stars[&tr.v[0]].pt,
                    &x,
                    self.robust_predicates,
                ) != -1
                {
                    break;
                } else {
                    //-- walk to incident to tr1,tr2
                    // println!("here");
                    let prev = self.stars[&tr.v[2]].link.get_prev_vertex(tr.v[0]).unwrap();
                    tr.v[1] = tr.v[2];
                    tr.v[2] = prev;
                }
            } else {
                //-- walk to incident to tr1,tr2
                // a.iter().position(|&x| x == 2), Some(1)
                let prev = self.stars[&tr.v[1]].link.get_prev_vertex(tr.v[2]).unwrap();
                tr.v[0] = tr.v[2];
                tr.v[2] = prev;
            }
        }
        return tr;
    }

    fn flip22(&mut self, tr: &Triangle, opposite: usize) -> (Triangle, Triangle) {
        //-- step 1.
        self.stars
            .get_mut(&tr.v[0])
            .unwrap()
            .link
            .insert_after_v(opposite, tr.v[1]);
        //-- step 2.
        self.stars.get_mut(&tr.v[1]).unwrap().link.delete(tr.v[2]);
        //-- step 3.
        self.stars
            .get_mut(&opposite)
            .unwrap()
            .link
            .insert_after_v(tr.v[0], tr.v[2]);
        //-- step 4.
        self.stars.get_mut(&tr.v[2]).unwrap().link.delete(tr.v[1]);
        //-- make 2 triangles to return (to stack)
        let ret0 = Triangle {
            v: [tr.v[0], tr.v[1], opposite],
        };
        let ret1 = Triangle {
            v: [tr.v[0], opposite, tr.v[2]],
        };
        (ret0, ret1)
    }

    fn get_opposite_vertex(&self, tr: &Triangle) -> usize {
        self.stars[&tr.v[2]].link.get_next_vertex(tr.v[1]).unwrap()
    }

    /// Returns a Vec<Vec<f64>> of all the vertices (including the infinite one)
    pub fn all_vertices(&self) -> Vec<Vec<f64>> {
        let mut pts: Vec<Vec<f64>> = Vec::with_capacity(self.stars.len() - 1);
        for i in 0..self.stars.len() {
            pts.push(self.stars.get(&i).unwrap().pt.to_vec());
        }
        pts
    }

    /// Returns a <Vec<usize> of all the finite edges (implicitly grouped by 2)
    pub fn all_edges(&self) -> Vec<usize> {
        let mut edges: Vec<usize> = Vec::new();
        for i in 1..self.stars.len() {
            for value in self.stars[&i].link.iter() {
                if (*value != 0) && (i < *value) {
                    edges.push(i);
                    edges.push(*value);
                }
            }
        }
        edges
    }

    /// Returns a <Vec<Triangle> of all the finite triangles (including the infinite one)
    pub fn all_triangles(&self) -> Vec<Triangle> {
        let mut trs: Vec<Triangle> = Vec::new();
        for (i, star) in &self.stars {
            //-- reconstruct triangles
            for (j, value) in star.link.iter().enumerate() {
                if i < value {
                    // let k = star.l[self.nexti(star.link.len(), j)];
                    let k = star.link[star.link.next_index(j)];
                    if i < &k {
                        let tr = Triangle { v: [*i, *value, k] };
                        if tr.is_infinite() == false {
                            // println!("{}", tr);
                            trs.push(tr);
                        }
                    }
                }
            }
        }
        trs
    }

    pub fn all_triangles_active(&self) -> Vec<Triangle> {
        let mut trs: Vec<Triangle> = Vec::new();
        for (i, star) in &self.stars {
            //-- reconstruct triangles
            for (j, value) in star.link.iter().enumerate() {
                if (i < value) && (self.stars.contains_key(&value) == true) {
                    // let k = star.l[self.nexti(star.link.len(), j)];
                    let k = star.link[star.link.next_index(j)];
                    if (i < &k) && (self.stars.contains_key(&k) == true) {
                        let tr = Triangle { v: [*i, *value, k] };
                        if tr.is_infinite() == false {
                            // println!("{}", tr);
                            trs.push(tr);
                        }
                    }
                }
            }
        }
        trs
    }

    /// Validates the Delaunay triangulation:
    /// (1) checks each triangle against each vertex (circumcircle tests); very slow
    /// (2) checks whether the convex hull is really convex
    pub fn is_valid(&self) -> bool {
        self.is_valid_ch_convex() && self.is_valid_circumcircle()
    }

    fn is_valid_circumcircle(&self) -> bool {
        let mut re = true;
        let trs = self.all_triangles();
        for tr in trs.iter() {
            for i in 1..self.stars.len() {
                if geom::incircle(
                    &self.stars[&tr.v[0]].pt,
                    &self.stars[&tr.v[1]].pt,
                    &self.stars[&tr.v[2]].pt,
                    &self.stars[&i].pt,
                    self.robust_predicates,
                ) > 0
                {
                    println!("NOT DELAUNAY FFS!");
                    println!("{} with {}", tr, i);
                    re = false
                }
            }
        }
        re
    }

    fn is_valid_ch_convex(&self) -> bool {
        let mut re = true;
        let ch = self.convex_hull();
        for i in 0..ch.len() {
            if geom::orient2d(
                &self.stars[&ch[i % ch.len()]].pt,
                &self.stars[&ch[(i + 1) % ch.len()]].pt,
                &self.stars[&ch[(i + 2) % ch.len()]].pt,
                self.robust_predicates,
            ) == -1
            {
                re = false;
                break;
            }
        }
        if re == false {
            println!("CONVEX NOT CONVEX");
        }
        return re;
    }

    pub fn printme(&self, withxyz: bool) -> String {
        let mut s = String::from("**********\n");
        // s.push_str(&format!("#pts: {}\n", self.number_pts()));

        let mut allkeys: Vec<&usize> = self.stars.keys().collect();
        allkeys.sort();
        for each in allkeys {
            let v = self.stars.get(each).unwrap();
            s.push_str(&format!("{}: [", *each));
            for each2 in v.link.iter() {
                s.push_str(&format!("{} - ", each2));
            }
            s.push_str(&format!("]\n"));
        }

        // for (i, p) in &self.stars {
        //     // for (i, p) in self.stars.iter().enumerate() {
        //     // s.push_str(&format!("{}: {}\n", i, self.stars[i].link));
        //     s.push_str(&format!("{}: [", i));
        //     for each in p.link.iter() {
        //         s.push_str(&format!("{} - ", each));
        //     }
        //     s.push_str(&format!("]\n"));
        //     if withxyz == true {
        //         s.push_str(&format!("\t{:?}\n", self.stars[&i].pt));
        //     }
        // }
        s.push_str("**********\n");
        s
    }

    fn vertex_exists(&self, v: usize) -> bool {
        self.stars.contains_key(&v)
    }

    /// write a GeoJSON file of the triangles/vertices to disk
    pub fn write_geojson_triangles(&self, path: String) -> std::io::Result<()> {
        let mut fc = FeatureCollection {
            bbox: None,
            features: vec![],
            foreign_members: None,
        };
        //-- vertices
        for (i, star) in &self.stars {
            if *i == 0 {
                continue;
            }
            let pt = Geometry::new(Value::Point(vec![star.pt[0], star.pt[1]]));
            let mut attributes = Map::new();
            attributes.insert(String::from("id"), to_value(i.to_string()).unwrap());
            attributes.insert(
                String::from("active"),
                serde_json::value::Value::Bool(star.active),
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
        for (i, star) in &self.stars {
            for (j, value) in star.link.iter().enumerate() {
                if (i < value) && (self.stars.contains_key(&value) == true) {
                    let k = star.link[star.link.next_index(j)];
                    if (i < &k) && (self.stars.contains_key(&k) == true) {
                        let tr = Triangle { v: [*i, *value, k] };
                        if tr.is_infinite() == false {
                            let mut l: Vec<Vec<Vec<f64>>> = vec![vec![Vec::with_capacity(1); 4]];
                            l[0][0].push(self.stars[i].pt[0]);
                            l[0][0].push(self.stars[i].pt[1]);
                            l[0][1].push(self.stars[value].pt[0]);
                            l[0][1].push(self.stars[value].pt[1]);
                            l[0][2].push(self.stars[&k].pt[0]);
                            l[0][2].push(self.stars[&k].pt[1]);
                            l[0][3].push(self.stars[i].pt[0]);
                            l[0][3].push(self.stars[i].pt[1]);
                            let gtr = Geometry::new(Value::Polygon(l));
                            // let mut attributes = Map::new();
                            // if self.stars[]
                            // attributes.insert(String::from("active"), to_value();

                            let f = Feature {
                                bbox: None,
                                geometry: Some(gtr),
                                id: None,
                                properties: None, //Some(attributes),
                                foreign_members: None,
                            };
                            fc.features.push(f);
                        }
                    }
                }
            }
        }

        //-- write the file to disk
        let mut fo = File::create(path)?;
        write!(fo, "{}", fc.to_string()).unwrap();
        Ok(())
    }

    /// write a GeoJSON file of the quadtree/grid to disk
    pub fn write_geojson_grid(&self, path: String) -> std::io::Result<()> {
        let mut fc = FeatureCollection {
            bbox: None,
            features: vec![],
            foreign_members: None,
        };
        for (i, w) in self.qt.gpts.iter().enumerate() {
            for (j, h) in w.iter().enumerate() {
                let mut l: Vec<Vec<Vec<f64>>> = vec![vec![Vec::with_capacity(1); 5]];
                l[0][0].push(self.qt.minx + ((i * self.qt.cellsize) as f64));
                l[0][0].push(self.qt.miny + ((j * self.qt.cellsize) as f64));
                l[0][1].push(self.qt.minx + (((i + 1) * self.qt.cellsize) as f64));
                l[0][1].push(self.qt.miny + ((j * self.qt.cellsize) as f64));
                l[0][2].push(self.qt.minx + (((i + 1) * self.qt.cellsize) as f64));
                l[0][2].push(self.qt.miny + (((j + 1) * self.qt.cellsize) as f64));
                l[0][3].push(self.qt.minx + ((i * self.qt.cellsize) as f64));
                l[0][3].push(self.qt.miny + (((j + 1) * self.qt.cellsize) as f64));
                l[0][4].push(self.qt.minx + ((i * self.qt.cellsize) as f64));
                l[0][4].push(self.qt.miny + ((j * self.qt.cellsize) as f64));

                let gtr = Geometry::new(Value::Polygon(l));
                let f = Feature {
                    bbox: None,
                    geometry: Some(gtr),
                    id: None,
                    properties: None, //Some(attributes),
                    foreign_members: None,
                };
                fc.features.push(f);
            }
        }
        //-- write the file to disk
        let mut fo = File::create(path)?;
        write!(fo, "{}", fc.to_string()).unwrap();
        Ok(())
    }
}

// for (i, star) in &self.stars {
//     //-- reconstruct triangles
//     for (j, value) in star.link.iter().enumerate() {
//         if (i < value) && (self.stars.contains_key(&value) == true) {
//             // let k = star.l[self.nexti(star.link.len(), j)];
//             let k = star.link[star.link.next_index(j)];
//             if (i < &k) && (self.stars.contains_key(&k) == true) {
//                 let tr = Triangle { v: [*i, *value, k] };
//                 if tr.is_infinite() == false {
//                     // println!("{}", tr);
//                     trs.push(tr);
//                 }
//             }
//         }
//     }
// }

// let trs = self.all_triangles_active();
// let mut f = File::create(path)?;
// let mut s = String::new();
// let mut ids: Vec<&usize> = self.stars.keys().collect();
// ids.sort();
// println!("{:?}", ids);

// for i in 1..self.stars.len() {
//     if twod == true {
//         s.push_str(&format!(
//             "v {} {} {}\n",
//             self.stars[&i].pt[0], self.stars[&i].pt[1], 0
//         ));
//     } else {
//         s.push_str(&format!(
//             "v {} {} {}\n",
//             self.stars[&i].pt[0], self.stars[&i].pt[1], self.stars[&i].pt[2]
//         ));
//     }
// }
// write!(f, "{}", s).unwrap();
// let mut s = String::new();
// for tr in trs.iter() {
//     s.push_str(&format!("f {} {} {}\n", tr.v[0], tr.v[1], tr.v[2]));
// }
// write!(f, "{}", s).unwrap();
// // println!("write fobj: {:.2?}", starttime.elapsed());
//         Ok(())
//     }
// }

impl fmt::Display for Triangulation {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        fmt.write_str("======== TRIANGULATION ========\n")?;
        fmt.write_str(&format!("# vertices: {:19}\n", self.number_of_vertices()))?;
        fmt.write_str(&format!("# triangles: {:18}\n", self.number_of_triangles()))?;
        fmt.write_str(&format!(
            "# convex hull: {:16}\n",
            self.number_of_vertices_on_convex_hull()
        ))?;
        fmt.write_str(&format!("---\nrobust: {}\n", self.robust_predicates))?;
        fmt.write_str("===============================\n")?;
        Ok(())
    }
}
