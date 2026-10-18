//-- ttview

//! A minimal 2D viewer for the `.sma` stream produced by `ttdt`.
//!
//! This is the standalone viewer (for a saved file or a remote pipe). For a
//! live view of `ttdt`'s internal state (quadtree cells + current
//! triangulation), use `ttdt --gui` instead.
//!
//! The stream is shown as a growing overview: one small dot per triangle, at
//! its centroid. This deliberately avoids drawing the triangles themselves — a
//! large LAZ file yields millions of triangles, and egui/wgpu pack all geometry
//! into a single vertex buffer, so only a bounded number of vertices can be
//! drawn at once. Dots are placed on a screen-space grid (at most one per
//! cell), which keeps the count bounded by the viewport and stable as the
//! stream grows.

use std::collections::HashSet;
use std::sync::mpsc::{Receiver, TryRecvError};
use std::sync::Arc;

use clap::Parser;
use eframe::egui;

use sst::view::{spawn_lines, Cam};

/// At most this many lines are applied to the scene per frame, so that a fast
/// producer cannot starve rendering.
const MAX_LINES_PER_FRAME: usize = 1_000_000;
/// Size of a dot, and of the screen-space grid cell holding (at most) one dot.
const DOT_PX: f32 = 2.0;
/// Flat colour of the dots.
const DOT_COLOR: egui::Color32 = egui::Color32::from_rgb(150, 180, 220);

#[derive(Parser)]
#[command(name = "ttview")]
#[command(about = "tintin -- 2D overview of the .sma stream of ttdt [ttview]")]
#[command(author, version)]
struct Cli {
    /// Optional .sma file to replay; reads stdin when omitted
    input: Option<String>,
    /// Records to release per ~16ms tick (0 = as fast as possible)
    #[arg(short, long, default_value_t = 0)]
    speed: usize,
}

/// One parsed record of the stream.
enum Event {
    Bbox([f64; 4]),
    Vertex(f64, f64),
    Triangle(u32, u32, u32),
    Finalise,
}

/// The geometry accumulated so far.
struct Scene {
    /// index 0 is a placeholder so that the 1-based ids in the stream map
    /// directly to indices
    verts: Vec<[f64; 2]>,
    tris: Vec<[u32; 3]>,
    bbox: Option<[f64; 4]>,
    n_finalised: usize,
    done: bool,
}

impl Scene {
    fn new() -> Scene {
        Scene {
            verts: vec![[0.0, 0.0]],
            tris: Vec::new(),
            bbox: None,
            n_finalised: 0,
            done: false,
        }
    }

    fn apply(&mut self, e: Event) {
        match e {
            Event::Bbox(b) => self.bbox = Some(b),
            Event::Vertex(x, y) => self.verts.push([x, y]),
            Event::Triangle(a, b, c) => {
                let n = self.verts.len() as u32;
                //-- ids are 1-based; index 0 is the placeholder, never a triangle
                if a >= 1 && b >= 1 && c >= 1 && a < n && b < n && c < n {
                    self.tris.push([a, b, c]);
                }
            }
            Event::Finalise => self.n_finalised += 1,
        }
    }
}

struct Viewer {
    scene: Scene,
    rx: Receiver<String>,

    cam: Cam,
    need_fit: bool,

    //-- the dot mesh, plus the state it was built from. Dots are only added for
    //-- screen cells that are not already occupied, so the count stays bounded
    //-- by the viewport; the mesh is extended incrementally as triangles arrive.
    mesh: Arc<egui::Mesh>,
    occupancy: HashSet<(i32, i32)>,
    built_tris: usize,
    built_cam: Cam,
    drawn_dots: usize,
}

impl Viewer {
    fn new(cc: &eframe::CreationContext<'_>, cli: Cli) -> Viewer {
        let rx = spawn_lines(cli.input, cli.speed, cc.egui_ctx.clone());
        Viewer {
            scene: Scene::new(),
            rx,
            cam: Cam::new(),
            need_fit: true,
            mesh: Arc::new(egui::Mesh::default()),
            occupancy: HashSet::new(),
            built_tris: 0,
            built_cam: Cam::new(),
            drawn_dots: 0,
        }
    }

    fn drain(&mut self) {
        let mut n = 0;
        loop {
            match self.rx.try_recv() {
                Ok(l) => {
                    if let Some(e) = parse_line(&l) {
                        self.scene.apply(e);
                    }
                    n += 1;
                    if n >= MAX_LINES_PER_FRAME {
                        break;
                    }
                }
                Err(TryRecvError::Empty) => break,
                Err(TryRecvError::Disconnected) => {
                    self.scene.done = true;
                    break;
                }
            }
        }
    }

    /// Add a dot at the centroid of every triangle not yet integrated, skipping
    /// screen cells that already have one (and anything off-screen).
    fn incorporate(&mut self, cam: Cam) {
        let from = self.built_tris;
        let total = self.scene.tris.len();
        if from >= total {
            return;
        }
        let verts = &self.scene.verts;
        let mesh = Arc::make_mut(&mut self.mesh);
        let half = DOT_PX * 0.5;
        for tri in &self.scene.tris[from..total] {
            let a = &verts[tri[0] as usize];
            let b = &verts[tri[1] as usize];
            let c = &verts[tri[2] as usize];
            let p = cam.to_screen((a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0);
            if !cam.rect.contains(p) {
                continue;
            }
            let cell = ((p.x / DOT_PX).floor() as i32, (p.y / DOT_PX).floor() as i32);
            if self.occupancy.insert(cell) {
                let base = mesh.vertices.len() as u32;
                mesh.colored_vertex(egui::Pos2::new(p.x - half, p.y - half), DOT_COLOR);
                mesh.colored_vertex(egui::Pos2::new(p.x + half, p.y - half), DOT_COLOR);
                mesh.colored_vertex(egui::Pos2::new(p.x + half, p.y + half), DOT_COLOR);
                mesh.colored_vertex(egui::Pos2::new(p.x - half, p.y + half), DOT_COLOR);
                mesh.add_triangle(base, base + 1, base + 2);
                mesh.add_triangle(base, base + 2, base + 3);
                self.drawn_dots += 1;
            }
        }
        self.built_tris = total;
    }
}

impl eframe::App for Viewer {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        self.drain();

        ui.horizontal(|ui| {
            if ui.button("Fit view").clicked() {
                self.need_fit = true;
            }
            ui.label(format!(
                "vertices: {}",
                self.scene.verts.len().saturating_sub(1)
            ));
            ui.label(format!("triangles: {}", self.scene.tris.len()));
            ui.label(format!("finalised: {}", self.scene.n_finalised));
            if self.scene.done {
                ui.label("finished");
            } else {
                ui.label("streaming...");
            }
        });
        ui.separator();

        let (response, painter) =
            ui.allocate_painter(ui.available_size(), egui::Sense::click_and_drag());
        let rect = response.rect;

        let mut cam = self.cam;
        cam.rect = rect;

        //-- fit the whole bbox into the viewport once it is known
        if self.need_fit {
            if let Some(b) = self.scene.bbox {
                cam.fit_bbox(b);
                self.need_fit = false;
            }
        }
        cam.interact(ui, &response);
        self.cam = cam;

        //-- background
        painter.rect_filled(rect, 0.0, egui::Color32::from_gray(20));

        //-- if the camera moved, the dots must be placed from scratch
        if cam != self.built_cam {
            self.mesh = Arc::new(egui::Mesh::default());
            self.occupancy.clear();
            self.built_tris = 0;
            self.drawn_dots = 0;
            self.built_cam = cam;
        }

        //-- integrate any newly received triangles
        self.incorporate(cam);

        if self.mesh.indices.is_empty() == false {
            painter.add(egui::Shape::mesh(self.mesh.clone()));
        }

        //-- keep animating while the stream is still being read
        if !self.scene.done {
            ui.ctx().request_repaint();
        }
    }
}

/// Parse one line of the stream; `None` for comments, blank lines, and records
/// that the viewer does not need (`n`, `c`, `s`, `w`).
fn parse_line(l: &str) -> Option<Event> {
    let mut it = l.split_whitespace();
    match it.next()? {
        "#" => None,
        "b" => {
            let a = it.next()?.parse::<f64>().ok()?;
            let b = it.next()?.parse::<f64>().ok()?;
            let c = it.next()?.parse::<f64>().ok()?;
            let d = it.next()?.parse::<f64>().ok()?;
            Some(Event::Bbox([a, b, c, d]))
        }
        "v" => {
            let x = it.next()?.parse::<f64>().ok()?;
            let y = it.next()?.parse::<f64>().ok()?;
            Some(Event::Vertex(x, y))
        }
        "f" => {
            let a = it.next()?.parse::<u32>().ok()?;
            let b = it.next()?.parse::<u32>().ok()?;
            let c = it.next()?.parse::<u32>().ok()?;
            Some(Event::Triangle(a, b, c))
        }
        "x" => {
            //-- a finalised vertex; only the count matters for the overview
            Some(Event::Finalise)
        }
        _ => None,
    }
}

fn main() -> eframe::Result {
    let cli = Cli::parse();
    let native_options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1200.0, 800.0])
            .with_title("ttview"),
        ..Default::default()
    };
    eframe::run_native(
        "ttview",
        native_options,
        Box::new(move |cc| Ok(Box::new(Viewer::new(cc, cli)))),
    )
}
