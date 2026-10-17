//-- ttview

//! A minimal 2D viewer for the `.sma` stream produced by `ttdt`.
//!
//! The stream is read on a background thread and shown as a growing overview:
//! one small dot per triangle, at its centroid. This deliberately avoids
//! drawing the triangles themselves — a large LAZ file yields millions of
//! triangles, and egui/wgpu pack all geometry into a single vertex buffer, so
//! only a bounded number of vertices can be drawn at once. Dots are placed on a
//! screen-space grid (at most one per cell), which keeps the count bounded by
//! the viewport and stable as the stream grows.
//!
//! Reads from stdin by default (pipe `ttdt` into it), or replays a saved `.sma`
//! file. `--speed` paces the release of records so the streaming order can be
//! watched offline.

use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::sync::mpsc::{self, Receiver, TryRecvError};
use std::sync::Arc;
use std::thread;
use std::time::Duration;

use clap::Parser;
use eframe::egui;

/// With `--speed 0`, still repaint every this many records so that a very long
/// stream stays visible while it is being drained as fast as possible.
const REPAINT_EVERY: usize = 500_000;
/// Delay between batches when `--speed` is set.
const TICK: Duration = Duration::from_millis(16);
/// At most this many events are applied to the scene per frame, so that a fast
/// producer cannot starve rendering.
const MAX_EVENTS_PER_FRAME: usize = 1_000_000;
/// Bound on the reader->UI channel; provides backpressure so a fast reader
/// cannot buffer the whole stream in memory.
const CHANNEL_CAP: usize = 500_000;

/// Size of a dot, and of the screen-space grid cell that holds (at most) one
/// dot, in points.
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
    Done,
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
                if a < n && b < n && c < n {
                    self.tris.push([a, b, c]);
                }
            }
            Event::Finalise => self.n_finalised += 1,
            Event::Done => self.done = true,
        }
    }
}

/// The camera: the world point shown at the viewport centre, the number of
/// pixels per world unit, and a small screen-space pan offset.
#[derive(Clone, Copy, PartialEq)]
struct Cam {
    rect: egui::Rect,
    center: [f64; 2],
    scale: f32,
    pan: egui::Vec2,
}

impl Cam {
    fn to_screen(&self, x: f64, y: f64) -> egui::Pos2 {
        egui::Pos2::new(
            self.rect.center().x + self.pan.x + ((x - self.center[0]) as f32) * self.scale,
            self.rect.center().y + self.pan.y - ((y - self.center[1]) as f32) * self.scale,
        )
    }

    fn to_world(&self, p: egui::Pos2) -> [f64; 2] {
        [
            self.center[0] + ((p.x - self.rect.center().x - self.pan.x) / self.scale) as f64,
            self.center[1] - ((p.y - self.rect.center().y - self.pan.y) / self.scale) as f64,
        ]
    }
}

struct Viewer {
    scene: Scene,
    rx: Receiver<Event>,

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
        let rx = spawn_reader(cli.input, cli.speed, cc.egui_ctx.clone());
        let cam = Cam {
            rect: egui::Rect::NOTHING,
            center: [0.0, 0.0],
            scale: 1.0,
            pan: egui::Vec2::ZERO,
        };
        Viewer {
            scene: Scene::new(),
            rx,
            cam,
            need_fit: true,
            mesh: Arc::new(egui::Mesh::default()),
            occupancy: HashSet::new(),
            built_tris: 0,
            built_cam: cam,
            drawn_dots: 0,
        }
    }

    fn drain(&mut self) {
        let mut n = 0;
        loop {
            match self.rx.try_recv() {
                Ok(e) => {
                    self.scene.apply(e);
                    n += 1;
                    if n >= MAX_EVENTS_PER_FRAME {
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
                mesh.colored_vertex(egui::pos2(p.x - half, p.y - half), DOT_COLOR);
                mesh.colored_vertex(egui::pos2(p.x + half, p.y - half), DOT_COLOR);
                mesh.colored_vertex(egui::pos2(p.x + half, p.y + half), DOT_COLOR);
                mesh.colored_vertex(egui::pos2(p.x - half, p.y + half), DOT_COLOR);
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
                let w = b[2] - b[0];
                let h = b[3] - b[1];
                if w > 0.0 && h > 0.0 {
                    cam.scale = (rect.width() as f64 * 0.95 / w)
                        .min(rect.height() as f64 * 0.95 / h) as f32;
                    cam.center = [(b[0] + b[2]) / 2.0, (b[1] + b[3]) / 2.0];
                    cam.pan = egui::Vec2::ZERO;
                    self.need_fit = false;
                }
            }
        }

        //-- pan by dragging
        if response.dragged() {
            cam.pan += response.drag_delta();
        }

        //-- zoom about the cursor with the wheel
        if response.hovered() {
            let scroll = ui.input(|i| i.smooth_scroll_delta.y);
            if scroll != 0.0 {
                if let Some(pos) = response.hover_pos() {
                    let before = cam.to_world(pos);
                    cam.scale = (cam.scale * (scroll * 0.001).exp()).clamp(1e-9, 1e12);
                    cam.center[0] =
                        before[0] - ((pos.x - rect.center().x - cam.pan.x) as f64 / cam.scale as f64);
                    cam.center[1] =
                        before[1] + ((pos.y - rect.center().y - cam.pan.y) as f64 / cam.scale as f64);
                }
            }
        }
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

        painter.add(egui::Shape::mesh(self.mesh.clone()));

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

/// Read the stream on a background thread, sending parsed events to the UI.
fn spawn_reader(input: Option<String>, speed: usize, ctx: egui::Context) -> Receiver<Event> {
    let (tx, rx) = mpsc::sync_channel(CHANNEL_CAP);
    thread::spawn(move || {
        let source: Box<dyn BufRead> = match input {
            Some(path) => match File::open(&path) {
                Ok(f) => Box::new(BufReader::new(f)),
                Err(e) => {
                    eprintln!("ttview: cannot open {path}: {e}");
                    let _ = tx.send(Event::Done);
                    ctx.request_repaint();
                    return;
                }
            },
            None => Box::new(BufReader::new(io::stdin())),
        };

        let mut in_batch: usize = 0;
        for line in source.lines() {
            let l = match line {
                Ok(l) => l,
                Err(_) => break,
            };
            if let Some(ev) = parse_line(&l) {
                if tx.send(ev).is_err() {
                    //-- receiver gone: the window was closed
                    return;
                }
            }
            in_batch += 1;
            let flush = if speed > 0 {
                in_batch >= speed
            } else {
                in_batch >= REPAINT_EVERY
            };
            if flush {
                in_batch = 0;
                ctx.request_repaint();
                if speed > 0 {
                    thread::sleep(TICK);
                }
            }
        }
        let _ = tx.send(Event::Done);
        ctx.request_repaint();
    });
    rx
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
