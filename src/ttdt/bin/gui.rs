//-- ttdt --gui: a live 2D overview of the in-memory triangulation.

use std::sync::mpsc::{Receiver, TryRecvError};
use std::sync::Arc;
use std::time::{Duration, Instant};

use eframe::egui;
use eframe::egui::{Color32, Shape};

use sst::view::{push_line, spawn_lines, Cam};

use crate::handle_line;
use crate::triangulator::Triangulation;

/// Per-frame budget for feeding input to the triangulator, so that a fast
/// producer cannot starve rendering.
const BUDGET: Duration = Duration::from_millis(12);
/// Non-finalised cells (light blue, translucent).
const ACTIVE_CELL: Color32 = Color32::from_rgba_unmultiplied_const(120, 175, 240, 110);
/// Finalised cells (white; blends into the background).
const FINAL_CELL: Color32 = Color32::WHITE;
/// Background.
const BG: Color32 = Color32::WHITE;
/// Current triangulation (grey wireframe).
const WIRE: Color32 = Color32::from_gray(90);
const WIRE_W: f32 = 1.0;

pub fn run() -> std::io::Result<()> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1200.0, 800.0])
            .with_title("ttdt"),
        ..Default::default()
    };
    eframe::run_native("ttdt", options, Box::new(|cc| Ok(Box::new(View::new(cc)))))
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))
}

struct View {
    dt: Triangulation,
    rx: Receiver<String>,

    cam: Cam,
    need_fit: bool,

    //-- cached cell mesh, rebuilt only when the partition or camera changes
    cells: Arc<egui::Mesh>,
    built_rev: u64,
    built_cam: Cam,

    done: bool,
}

impl View {
    fn new(cc: &eframe::CreationContext<'_>) -> View {
        let rx = spawn_lines(None, 0, cc.egui_ctx.clone());
        View {
            dt: Triangulation::new(),
            rx,
            cam: Cam::new(),
            need_fit: true,
            cells: Arc::new(egui::Mesh::default()),
            built_rev: u64::MAX,
            built_cam: Cam::new(),
            done: false,
        }
    }

    /// Feed up to `BUDGET` worth of input to the triangulator.
    fn feed(&mut self) {
        let start = Instant::now();
        loop {
            match self.rx.try_recv() {
                Ok(l) => {
                    let _ = handle_line(&mut self.dt, &l);
                }
                Err(TryRecvError::Empty) => break,
                Err(TryRecvError::Disconnected) => {
                    if !self.done {
                        let _ = self.dt.flush_output();
                        self.done = true;
                    }
                    break;
                }
            }
            if start.elapsed() >= BUDGET {
                break;
            }
        }
    }

    /// The quadtree cell partition: light blue for non-finalised cells, white
    /// for finalised ones.
    fn build_cells(&self, cam: Cam, rect: egui::Rect) -> egui::Mesh {
        let mut mesh = egui::Mesh::default();
        for c in self.dt.cell_views() {
            let r = egui::Rect::from_two_pos(
                cam.to_screen(c.bbox[0], c.bbox[1]),
                cam.to_screen(c.bbox[2], c.bbox[3]),
            );
            if !rect.intersects(r) {
                continue;
            }
            let color = if c.finalised { FINAL_CELL } else { ACTIVE_CELL };
            let base = mesh.vertices.len() as u32;
            mesh.colored_vertex(r.left_top(), color);
            mesh.colored_vertex(r.right_top(), color);
            mesh.colored_vertex(r.right_bottom(), color);
            mesh.colored_vertex(r.left_bottom(), color);
            mesh.add_triangle(base, base + 1, base + 2);
            mesh.add_triangle(base, base + 2, base + 3);
        }
        mesh
    }

    /// The current (non-finalised) triangulation, as a grey wireframe.
    fn build_wire(&self, cam: Cam, rect: egui::Rect) -> egui::Mesh {
        let mut mesh = egui::Mesh::default();
        let Some(bbox) = self.dt.bbox() else {
            return mesh;
        };
        for t in self.dt.active_triangles() {
            //-- the augmented triangulation uses an "infinite" vertex at
            //-- (-99999.9, -99999.9); it is normally already excluded, but skip
            //-- anything that falls outside the data extent so no such geometry
            //-- can ever be drawn
            if t.iter().any(|v| {
                v[0] < bbox[0] - 1.0
                    || v[0] > bbox[2] + 1.0
                    || v[1] < bbox[1] - 1.0
                    || v[1] > bbox[3] + 1.0
            }) {
                continue;
            }
            let p0 = cam.to_screen(t[0][0], t[0][1]);
            let p1 = cam.to_screen(t[1][0], t[1][1]);
            let p2 = cam.to_screen(t[2][0], t[2][1]);
            if !rect.intersects(egui::Rect::from_points(&[p0, p1, p2])) {
                continue;
            }
            push_line(&mut mesh, p0, p1, WIRE_W, WIRE);
            push_line(&mut mesh, p1, p2, WIRE_W, WIRE);
            push_line(&mut mesh, p2, p0, WIRE_W, WIRE);
        }
        mesh
    }
}

impl eframe::App for View {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        self.feed();

        ui.horizontal(|ui| {
            if ui.button("Fit view").clicked() {
                self.need_fit = true;
            }
            ui.label(format!("cells: {}", self.dt.number_of_cells()));
            ui.label(format!(
                "active triangles: {}",
                self.dt.number_of_active_triangles()
            ));
            ui.label(format!("vertices: {}", self.dt.number_of_vertices()));
            if self.done {
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
        if self.need_fit {
            if let Some(b) = self.dt.bbox() {
                cam.fit_bbox(b);
                self.need_fit = false;
            }
        }
        cam.interact(ui, &response);
        self.cam = cam;

        //-- background
        painter.rect_filled(rect, 0.0, BG);

        //-- nothing meaningful can be placed before the bbox is known
        if self.dt.bbox().is_some() {
            //-- quadtree cells (cached; rebuilt when partition/camera changes)
            if cam != self.built_cam || self.dt.qt_rev() != self.built_rev {
                self.cells = Arc::new(self.build_cells(cam, rect));
                self.built_cam = cam;
                self.built_rev = self.dt.qt_rev();
            }
            if self.cells.indices.is_empty() == false {
                painter.add(Shape::mesh(self.cells.clone()));
            }

            //-- current triangulation (grey wireframe)
            let wire = self.build_wire(cam, rect);
            if wire.indices.is_empty() == false {
                painter.add(Shape::mesh(Arc::new(wire)));
            }
        }

        //-- keep animating while the stream is still being read
        if !self.done {
            ui.ctx().request_repaint();
        }
    }
}
