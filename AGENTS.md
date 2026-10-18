# AGENTS.md

## What this is
Rust implementation (Cargo package `sst`, toolset "tintin") of a **streaming Delaunay triangulation pipeline for LiDAR point clouds**, inspired by Isenburg's streaming triangulator. Algorithm reference: `help/Isenburg06.pdf`.

**Trust `Cargo.toml` over `README.md`.** The README predates the `sst*` → `tt*` rename: it still names old binaries (`sstdt`, `sstats`), says "three binaries", and claims `ttfin` reads XYZ (it does not).

## Build
```bash
cargo build --release   # the only build that is used; performance-critical
```
Binaries land in `target/release/`. There are no tests, no CI, and no rustfmt/clippy config.

## Binaries
| bin | input | output |
|-----|-------|--------|
| `ttfin <file-or-.files> <cellsize>` | LAS/LAZ on disk | `.spa` stream → stdout |
| `ttdt [--gui]` | `.spa` on **stdin** | `.sma` (OBJ-like) → stdout (and/or GUI) |
| `ttsimpl <vepsilon>` | `.spa` on **stdin** | `.spa` → stdout (thinning filter) |
| `ttstats <file>` | LAS/LAZ/`.txt`/`.xyz` | per-cell stats → stdout |
| `ttval` | `.sma` on **stdin** | validation counts → stdout |
| `ttview [file]` | `.sma` on **stdin** (or a file) | GUI window (dot overview) |

```bash
./target/release/ttfin ./data/crop.laz 10 | ./target/release/ttdt > ./data/crop.sma
# optional simplification stage between the two:
./target/release/ttfin ./data/crop.laz 10 | ./target/release/ttsimpl 0.5 | ./target/release/ttdt > ./data/crop.sma
```

Gotchas:
- `ttfin` accepts only `.las`/`.laz`/a `.files` list; other extensions print `input format not accepted`. Multiple files must be listed in a `.files` file that contains **absolute paths** (see `data/crop_1_2.files`), so it is machine-specific.
- `ttfin -g/--ground` keeps only LAS points with `Classification::Ground` (value 2). It is applied in passes 2 and 3 while `count` still advances over *all* points, so the sprinkled indices stay aligned with the header point count (pass 1 stays header-only). `ttstats` has no such filter.
- `ttstats` parses `std::env::args` directly (no clap): it has no `--help` and no verbosity flag, unlike the others.
- `ttdt`, `ttsimpl`, `ttval` only read stdin; a file-path argument is ignored.
- `ttdt --gui` shows a live view of `ttdt`'s **internal** state (quadtree cells + the current, non-finalised triangulation) — it still writes `.sma` to stdout. The GUI runs on the main thread (eframe requirement); stdin is read on a background thread and fed in under a per-frame time budget. The triangulation accessors it uses (`cell_views`, `active_triangles`, `number_of_cells`, `qt_rev`, `bbox`) were added to `triangulator.rs`.
- `ttview` (`--speed <N>` paces replay) draws **one dot per triangle** on a screen-space grid, not filled triangles: egui-wgpu packs *all* meshes into a single vertex buffer with a 256 MiB cap, so drawing millions of triangles panics — the dot overview keeps the vertex count bounded by the viewport.
- `ttview` and the `ttdt --gui` view are the only GUI code; both need a graphical desktop and pull a large dependency tree (`eframe`/`egui`, default `wgpu` backend). Shared viewer plumbing lives in the `sst::view` lib module. Note `eframe` 0.37's `App` trait uses `fn ui(&mut self, ui, frame)`, not the older `fn update(...)`.

## Stream format (`.spa` in, `.sma` out)
One record per line; the first character selects the type:
- `#` comment · `n` point count · `c` grid dimension · `s` cell size · `b minx miny maxx maxy`
- `v x y z` point · `w x y z` "sprinkled" point (hoisted to the top of the stream)
- `x gx gy` finalised cell (`.spa`) · `x id` finalised/dead vertex (`.sma`)
- `f a b c` triangle, 1-based ids (`.sma` output only)

Writing a stream by hand or synthesising one: **`c` must appear before any `v`/`x`** — `ttdt` initialises its quadtree on `c` and `insert_one_vi` panics on a missing key otherwise. `c` must also be a power of two (`ttfin` derives it as `2^ceil(log2(...))`). `.sma` vertex/triangle ids are 1-based.

## Layout
- `src/<bin>/bin/main.rs` — one entrypoint per binary (single Cargo package, not a workspace).
- `src/lib.rs` + `src/view/mod.rs` — shared library: the pan/zoom `Cam` and the background line reader (`spawn_lines`), used by `ttview` and `ttdt --gui`.
- `src/ttdt/bin/main.rs` (CLI) + `src/ttdt/bin/gui.rs` (`--gui` app) — `handle_line` in `main.rs` is shared by both paths.
- `src/ttdt/bin/triangulator.rs` (+ `triangulator/geom.rs`) — the triangulator, used **only by `ttdt`**; it is a forked/trimmed copy of the external `startin` crate, so edits here do not affect the other binaries. Its `--gui` accessors live here too.
- `ttsimpl` instead depends on the published `startin = 0.8` (2D DT; the viewer-era `interpolate_tin_linear` is now the batch `startin::interpolation::interpolate(&TIN {}, ...)` API); `ttdt` does not use `startin`.
- LAS/LAZ I/O uses the `las` crate (`laz` feature); predicates use `robust`.
- `data/` and `Cargo.lock` are in `.gitignore` yet present/tracked; the large LAZ files are samples, not for editing.
- `misc_tools/` — throwaway Python helpers.

## Repo state
Branch `reeng` is mid-rename `sst*` → `tt*`: `src/tt*/` are untracked and the old `src/sst*/` files show as deleted in `git status`. Do not "fix" the deleted paths or re-add the `sst*` binaries.
