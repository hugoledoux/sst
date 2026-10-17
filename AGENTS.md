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
| `ttdt` | `.spa` on **stdin** | `.sma` (OBJ-like) → stdout |
| `ttsimpl <vepsilon>` | `.spa` on **stdin** | `.spa` → stdout (thinning filter) |
| `ttstats <file>` | LAS/LAZ/`.txt`/`.xyz` | per-cell stats → stdout |
| `ttval` | `.sma` on **stdin** | validation counts → stdout |
| `ttview [file]` | `.sma` on **stdin** (or a file) | GUI window (2D triangles) |

```bash
./target/release/ttfin ./data/crop.laz 10 | ./target/release/ttdt > ./data/crop.sma
# optional simplification stage between the two:
./target/release/ttfin ./data/crop.laz 10 | ./target/release/ttsimpl 0.5 | ./target/release/ttdt > ./data/crop.sma
```

Gotchas:
- `ttfin` accepts only `.las`/`.laz`/a `.files` list; other extensions print `input format not accepted`. Multiple files must be listed in a `.files` file that contains **absolute paths** (see `data/crop_1_2.files`), so it is machine-specific.
- `ttstats` parses `std::env::args` directly (no clap): it has no `--help` and no verbosity flag, unlike the others.
- `ttdt`, `ttsimpl`, `ttval` only read stdin; a file-path argument is ignored.
- `ttview` is the only GUI binary (`eframe`/`egui`, default `wgpu` backend); it needs a graphical desktop and pulls in a large dependency tree. `--speed <N>` paces reading for replay. Note `eframe` 0.37's `App` trait uses `fn ui(&mut self, ui, frame)`, not the older `fn update(...)`. It draws **one dot per triangle** on a screen-space grid (not filled triangles): egui-wgpu packs *all* meshes into a single vertex buffer with a 256 MiB cap, so drawing millions of triangles panics — the dot overview keeps the vertex count bounded by the viewport.

## Stream format (`.spa` in, `.sma` out)
One record per line; the first character selects the type:
- `#` comment · `n` point count · `c` grid dimension · `s` cell size · `b minx miny maxx maxy`
- `v x y z` point · `w x y z` "sprinkled" point (hoisted to the top of the stream)
- `x gx gy` finalised cell (`.spa`) · `x id` finalised/dead vertex (`.sma`)
- `f a b c` triangle, 1-based ids (`.sma` output only)

Writing a stream by hand or synthesising one: **`c` must appear before any `v`/`x`** — `ttdt` initialises its quadtree on `c` and `insert_one_vi` panics on a missing key otherwise. `c` must also be a power of two (`ttfin` derives it as `2^ceil(log2(...))`). `.sma` vertex/triangle ids are 1-based.

## Layout
- `src/<bin>/bin/main.rs` — one entrypoint per binary (single Cargo package, not a workspace).
- `src/ttdt/bin/triangulator.rs` (+ `triangulator/geom.rs`) — the triangulator, used **only by `ttdt`**; it is a forked/trimmed copy of the external `startin` crate, so edits here do not affect the other binaries.
- `ttsimpl` instead depends on the published `startin = 0.6.2` (2D DT with `interpolate_tin_linear`); `ttdt` does not.
- LAS/LAZ I/O uses the `las` crate (`laz` feature); predicates use `robust`.
- `data/` and `Cargo.lock` are in `.gitignore` yet present/tracked; the large LAZ files are samples, not for editing.
- `misc_tools/` — throwaway Python helpers.

## Repo state
Branch `reeng` is mid-rename `sst*` → `tt*`: `src/tt*/` are untracked and the old `src/sst*/` files show as deleted in `git status`. Do not "fix" the deleted paths or re-add the `sst*` binaries.
