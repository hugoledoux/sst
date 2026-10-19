# AGENTS.md

## What this is
Rust implementation (Cargo package `sst`, toolset "tintin") of a **streaming Delaunay triangulation pipeline for LiDAR point clouds**, inspired by Isenburg's streaming triangulator. Algorithm reference: `help/Isenburg06.pdf`.

**Trust `Cargo.toml` over `README.md`.** The README predates the `sst*` → `tt*` rename: it still names old binaries (`sstdt`, `sstats`), says "three binaries", and claims `ttfin` reads XYZ (it does not).

## Build
```bash
cargo build --release   # the only build that is used; performance-critical
```
Binaries land in `target/release/`. There are no tests, no CI, and no rustfmt/clippy config. `Cargo.lock` is gitignored and untracked: in a fresh checkout or `git worktree`, cargo re-resolves `geojson = "*"` to an API-incompatible version and the build fails — copy a live `Cargo.lock` in first.

## Binaries
| bin | input | output |
|-----|-------|--------|
| `ttfin <file-or-.files> <cellsize>` | LAS/LAZ on disk | `.spa` stream → stdout |
| `ttdt [--gui]` | `.spa` on **stdin** | `.sma` (OBJ-like) → stdout (and/or GUI) |
| `ttsimpl <vepsilon>` | `.spa` on **stdin** | `.spa` → stdout (thinning filter) |
| `ttstats <file>` | LAS/LAZ/`.txt`/`.xyz` | per-cell stats → stdout |
| `ttval` | `.sma` on **stdin** | validation counts → stdout |
| `spb2obj` | `.sma` on **stdin** | Wavefront OBJ text (`v`/`f` records only) → stdout |
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
- `ttdt`, `ttsimpl`, `ttval`, `spb2obj` only read stdin; a file-path argument is ignored.
- `ttdt --gui` shows a live view of `ttdt`'s **internal** state (quadtree cells + the current, non-finalised triangulation) — it still writes `.sma` to stdout. The GUI runs on the main thread (eframe requirement); stdin is read on a background thread and fed in under a per-frame time budget. The triangulation accessors it uses (`cell_views`, `active_triangles`, `number_of_cells`, `qt_rev`, `bbox`) were added to `triangulator.rs`.
- `ttview` (`--speed <N>` paces replay) draws **one dot per triangle** on a screen-space grid, not filled triangles: egui-wgpu packs *all* meshes into a single vertex buffer with a 256 MiB cap, so drawing millions of triangles panics — the dot overview keeps the vertex count bounded by the viewport.
- `ttview` and the `ttdt --gui` view are the only GUI code; both need a graphical desktop and pull a large dependency tree (`eframe`/`egui`, default `wgpu` backend). Shared viewer plumbing lives in the `sst::view` lib module. Note `eframe` 0.37's `App` trait uses `fn ui(&mut self, ui, frame)`, not the older `fn update(...)`.

## Stream format (`.spa` in, `.sma` out)
One record per line; the first character selects the type:
- `#` comment · `n` point count · `c` grid dimension · `s` cell size · `b minx miny maxx maxy`
- `v x y z` point · `w x y z` "sprinkled" point (hoisted to the top of the stream)
- `x gx gy` finalised cell (`.spa`) · `x id` finalised/dead vertex (`.sma`)
- `f a b c` triangle, 1-based ids (`.sma` output only)

Both formats also have a **binary variant, the default for `ttfin`/`ttdt` output** (`--ascii` opts back into text; every reader auto-detects via the magic `TTSPA1\n` / `TTSMA1\n`, so no reader flag exists). Records are the ASCII ones with fixed little-endian payloads: tag byte + `n` u64 · `c` u32 · `s` u64 · `b` 4×f64 · `v`/`w` 3×f64 · `x` 2×u32 (`.spa`); `b` 4×f64 · `v` 3×f64 · `f` 3×u32 · `x` u32 (`.sma`). The single codec lives in `src/stream.rs` (`sst::stream`); `record.to_line()` converts a decoded record back to its ASCII line, which is how the GUI viewers consume binary streams unchanged. `ttsimpl` mirrors the format of its input. Binary carries full f64 coordinates; `ttfin --ascii` quantises to the LAS scale digits.

Writing a stream by hand or synthesising one: **`c` must appear before any `v`/`x`** — `ttdt` initialises its quadtree on `c` and `insert_one_vi` panics on a missing key otherwise. `c` must also be a power of two (`ttfin` derives it as `2^ceil(log2(...))`). `.sma` vertex/triangle ids are 1-based.

## Layout
- `src/<bin>/bin/main.rs` — one entrypoint per binary (single Cargo package, not a workspace).
- `src/lib.rs` + `src/stream.rs` + `src/view/mod.rs` — shared library: the `.spa`/`.sma` codecs (`sst::stream`), the pan/zoom `Cam`, and the background line reader (`spawn_lines`), used by `ttview` and `ttdt --gui`.
- `src/ttdt/bin/main.rs` (CLI) + `src/ttdt/bin/gui.rs` (`--gui` app) — `handle_line` (ASCII) and `handle_record` (binary) in `main.rs` are shared by both paths.
- `src/ttdt/bin/triangulator.rs` (+ `triangulator/geom.rs`) — the triangulator, used **only by `ttdt`**; it is a forked/trimmed copy of the external `startin` crate, so edits here do not affect the other binaries. Its `--gui` accessors live here too.
- `ttsimpl` instead depends on the published `startin = 0.8` (2D DT; the viewer-era `interpolate_tin_linear` is now the batch `startin::interpolation::interpolate(&TIN {}, ...)` API); `ttdt` does not use `startin`.
- LAS/LAZ I/O uses the `las` crate (`laz` feature); predicates use `robust`.
- `data/` is gitignored, but the small samples (`crop.laz`, `crop_1.laz`, `crop_2.laz`, `crop_1_2.files`) are force-tracked; the `*.xyz` and large LAZ files exist only locally. Samples are inputs, not for editing.
- `misc_tools/` — throwaway Python helpers.

## Verifying changes
There are no tests; the de facto validation of a pipeline change is:
- run it and check Euler: `T = 2n − 2 − h` (h = convex-hull vertices, a small plausible number), plus a duplicate-face check;
- compare outputs **geometrically** — sorted coordinate triples per face — never by `.sma` vertex ids (they follow activation order) or byte-wise;
- `ttfin` re-randomises its sprinkles on every run, so two runs on the same input differ marginally (the 1 mm snap-tolerance lottery can shift `n` by a few points). To compare two binaries, feed them the **same** `.spa` file — save one with `ttfin > file.spa` first;
- quick input: `data/crop.laz` (~212k pts). Larger local-only samples: `AHN5_C_192000_309000.laz` (26M, ground-filtered; ~3.0 flips/insertion) and `37EN2_11.LAZ` (51.9M, raw era; ~3.5).

Verbosity (clap-verbosity-flag): default = errors only; `-v` = WARN — the `walk.3`/`walk.4` fallback warnings, whose counts are the health signal of the triangulator's walk; `-vv` = INFO — per-cell finalisation logs and the counters (`# insertions`, `# flips22`, `# walks`/`# walk steps`).

## Performance notes
- The pipeline's cost is robust predicates + finalisation + LAZ decompression; the streams are cheap *only because* every writer uses a `BufWriter` — per-line `write_all` once cost 208 s of system time on a 26M-point run. Never reintroduce them.
- `ttfin <cellsize>` is in metres; the grid becomes `2^ceil(log2(extent/cellsize))` per side (the `c` record). Aim for a few thousand points per cell (≈10 m at AHN5 density): larger cells make the per-cell finalisation batch exceed CPU caches (wall time grows superlinearly per point); smaller cells multiply chunk-boundary walk fallbacks (`walk.4` brute-force scans).
- Load-bearing invariant in `triangulator.rs`: `is_triangle_final` certifies only triangles whose 3 vertices lie in the finalised cell's bbox — that is what keeps the stars of *unfinalised*-cell vertices intact and alive. The walk fallbacks, star collection, and the `vs_dead` mark-and-sweep all rely on it.

## Repo state
The `sst*` → `tt*` rename is complete and committed; the only remaining `sst` is the Cargo package/lib name (binaries do `use sst::stream` / `sst::view`). `help/` (the reference paper) is untracked. The README still uses the old binary names — see the note at the top.
