//! Shared plumbing for the 2D viewers (`ttview` and the `ttdt --gui` view):
//! a simple pan/zoom camera and a background line reader.

use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::sync::mpsc::{self, Receiver};
use std::thread;
use std::time::Duration;

use eframe::egui;

/// Repaint every this many lines when reading as fast as possible.
const REPAINT_EVERY: usize = 500_000;
/// Delay between batches when a speed is set.
const TICK: Duration = Duration::from_millis(16);
/// Bound on the reader->UI channel; provides backpressure so a fast reader
/// cannot buffer the whole stream in memory.
const CHANNEL_CAP: usize = 500_000;

/// A 2D camera: the world point shown at the viewport centre, the number of
/// pixels per world unit, and a small screen-space pan offset.
#[derive(Clone, Copy, PartialEq)]
pub struct Cam {
    pub rect: egui::Rect,
    pub center: [f64; 2],
    pub scale: f32,
    pub pan: egui::Vec2,
}

impl Cam {
    pub fn new() -> Cam {
        Cam {
            rect: egui::Rect::NOTHING,
            center: [0.0, 0.0],
            scale: 1.0,
            pan: egui::Vec2::ZERO,
        }
    }

    pub fn to_screen(&self, x: f64, y: f64) -> egui::Pos2 {
        egui::Pos2::new(
            self.rect.center().x + self.pan.x + ((x - self.center[0]) as f32) * self.scale,
            self.rect.center().y + self.pan.y - ((y - self.center[1]) as f32) * self.scale,
        )
    }

    pub fn to_world(&self, p: egui::Pos2) -> [f64; 2] {
        [
            self.center[0] + ((p.x - self.rect.center().x - self.pan.x) / self.scale) as f64,
            self.center[1] - ((p.y - self.rect.center().y - self.pan.y) / self.scale) as f64,
        ]
    }

    /// Fit a world bbox to the viewport (95% margin).
    pub fn fit_bbox(&mut self, bbox: [f64; 4]) {
        let w = bbox[2] - bbox[0];
        let h = bbox[3] - bbox[1];
        if w > 0.0 && h > 0.0 {
            self.scale = (self.rect.width() as f64 * 0.95 / w)
                .min(self.rect.height() as f64 * 0.95 / h) as f32;
            self.center = [(bbox[0] + bbox[2]) / 2.0, (bbox[1] + bbox[3]) / 2.0];
            self.pan = egui::Vec2::ZERO;
        }
    }

    /// Pan by dragging and zoom about the cursor with the wheel.
    pub fn interact(&mut self, ui: &egui::Ui, response: &egui::Response) {
        if response.dragged() {
            self.pan += response.drag_delta();
        }
        if response.hovered() {
            let scroll = ui.input(|i| i.smooth_scroll_delta.y);
            if scroll != 0.0 {
                if let Some(pos) = response.hover_pos() {
                    let before = self.to_world(pos);
                    self.scale = (self.scale * (scroll * 0.001).exp()).clamp(1e-9, 1e12);
                    self.center[0] = before[0]
                        - ((pos.x - self.rect.center().x - self.pan.x) as f64 / self.scale as f64);
                    self.center[1] = before[1]
                        + ((pos.y - self.rect.center().y - self.pan.y) as f64 / self.scale as f64);
                }
            }
        }
    }
}

/// Append a thick line segment to an `egui::Mesh` (as a quad).
pub fn push_line(
    mesh: &mut egui::Mesh,
    a: egui::Pos2,
    b: egui::Pos2,
    thickness: f32,
    color: egui::Color32,
) {
    let d = b - a;
    let len = d.length();
    if len < 1e-6 {
        return;
    }
    let n = egui::Vec2::new(-d.y, d.x) / len * (thickness * 0.5);
    let base = mesh.vertices.len() as u32;
    mesh.colored_vertex(a + n, color);
    mesh.colored_vertex(b + n, color);
    mesh.colored_vertex(b - n, color);
    mesh.colored_vertex(a - n, color);
    mesh.add_triangle(base, base + 1, base + 2);
    mesh.add_triangle(base, base + 2, base + 3);
}

/// Read lines from `input` (stdin when `None`) on a background thread and send
/// them to the UI. `speed` is the number of lines released per ~16ms tick
/// (`0` = as fast as possible). Dropping the receiver stops the thread.
///
/// A binary stream (`.spa` or `.sma`) is decoded here and its records are sent
/// as the equivalent ASCII lines, so that the viewers work with either
/// format without further changes.
pub fn spawn_lines(input: Option<String>, speed: usize, ctx: egui::Context) -> Receiver<String> {
    let (tx, rx) = mpsc::sync_channel(CHANNEL_CAP);
    thread::spawn(move || {
        let mut source: Box<dyn BufRead> = match input {
            Some(path) => match File::open(&path) {
                Ok(f) => Box::new(BufReader::new(f)),
                Err(e) => {
                    eprintln!("cannot open {path}: {e}");
                    ctx.request_repaint();
                    return;
                }
            },
            None => Box::new(BufReader::new(io::stdin())),
        };

        let binary_spa = crate::stream::peek_is_binary(&mut source, crate::stream::SPA_MAGIC)
            .unwrap_or(false);
        let binary_sma = if binary_spa {
            false
        } else {
            crate::stream::peek_is_binary(&mut source, crate::stream::SMA_MAGIC).unwrap_or(false)
        };

        let mut in_batch: usize = 0;
        let mut tick = |in_batch: &mut usize| {
            let flush = if speed > 0 {
                *in_batch >= speed
            } else {
                *in_batch >= REPAINT_EVERY
            };
            if flush {
                *in_batch = 0;
                ctx.request_repaint();
                if speed > 0 {
                    thread::sleep(TICK);
                }
            }
        };

        if binary_spa {
            let _ = crate::stream::skip_magic(&mut source, crate::stream::SPA_MAGIC);
            while let Ok(Some(rec)) = crate::stream::spa_read(&mut source) {
                if tx.send(rec.to_line()).is_err() {
                    //-- receiver gone: the window was closed
                    return;
                }
                in_batch += 1;
                tick(&mut in_batch);
            }
        } else if binary_sma {
            let _ = crate::stream::skip_magic(&mut source, crate::stream::SMA_MAGIC);
            while let Ok(Some(rec)) = crate::stream::sma_read(&mut source) {
                if tx.send(rec.to_line()).is_err() {
                    //-- receiver gone: the window was closed
                    return;
                }
                in_batch += 1;
                tick(&mut in_batch);
            }
        } else {
            for line in source.lines() {
                let l = match line {
                    Ok(l) => l,
                    Err(_) => break,
                };
                if tx.send(l).is_err() {
                    //-- receiver gone: the window was closed
                    return;
                }
                in_batch += 1;
                tick(&mut in_batch);
            }
        }
        ctx.request_repaint();
    });
    rx
}
