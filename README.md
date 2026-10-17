# tintin

A Rust reimplementation (with new operators) of the streaming Delaunay triangulator of Martin Isenburg.

The streaming pipeline is (each stage reads/writes an ASCII stream, except the first):

```
LAS/LAZ files --> ttfin --> [ttsimpl] --> ttdt --> .sma (OBJ-like triangulation)
```

There are five binaries, all built together.


## To compile them all

```bash
cargo build --release
```

`--release` makes the whole pipeline a whoooole lot faster, but you can't debug.
The binaries end up in `./target/release/`.


## 1. **ttfin**

  - conceptually equivalent to Isenburg's *spfinalise*
  - finalises (chunks) a set of points in LAS/LAZ
  - can read several LAZ files at the same time, listed in a `.files` text file
    (absolute paths, one per line — see `./data/crop_1_2.files`)
  - reads the files on disk; writes the `.spa` stream to stdout

Usage: `ttfin <input.laz|input.files> <cellsize> [-g] [-s <sprinkle>]`

  - `<cellsize>`: the cell resolution as an integer, eg `10` for 10m×10m
  - `-s, --sprinkle`: fraction of points "sprinkled" (hoisted to the top of the
    stream), default `0.001`
  - `-g, --ground`: intended to use only ground points from the LAS file
    (accepted but not implemented yet)
  - verbosity is controlled with `-v`/`-vv`/`-q` (the INFO/WARN logs)

To see the format of the output stream ("spa" is **s**treaming **p**oint **a**scii; binary to follow at some point)

```bash
./target/release/ttfin -vv ./data/crop.laz 10 > ./data/crop_10.spa
```

Multiple LAZ files if you provide them in a `myfiles.files`, see the example in `./data/crop_1_2.files`.
This reads `crop_1.laz` and then `crop_2.laz` (in that order); I split `crop.laz` to obtain those 2 files.

```bash
./target/release/ttfin -vv ./data/crop_1_2.files 10 > /dev/null
```

(remove `-vv` to remove the INFO and WARN)

__.spa file__

```
# ttfin
n 212550
c 16
s 10
b 84600.000 447000.000 84699.998 447099.999
x 0 10
x 0 11
v 84601.790 447091.302 6.554
v 84660.438 447055.597 0.796
v 84637.450 447010.511 8.879
v 84662.152 447009.066 3.018
v 84669.030 447052.681 0.930
v 84603.720 447072.709 0.105
...
x 2 3
...
```

The records are: `n` (number of points), `c` (square grid dimension, a power of two), `s` (cell size), `b` (bbox), `x` (finalise a cell), `v` (vertex), `w` (sprinkled vertex).


## 2. **ttdt**

  - conceptually equivalent to Isenburg's *spdelaunay*
  - takes as input (stdin) the .spa output of ttfin and creates a DT of the points
  - outputs a variation of an .obj file where a few vertex finalisers (`x`) are added. It's simpler than Isenburg's .sma (**s**treaming **m**esh **a**scii) and can be opened as an .obj file by different software (MeshLab just ignores the extra lines, so you can open the file with it)

```bash
./target/release/ttfin -vv ./data/crop.laz 10 | ./target/release/ttdt -vv > ./data/crop_10.sma
```

The output file `crop_10.sma` is a valid OBJ file: change its extension and you can open it normally in [MeshLab](https://www.meshlab.net/) for instance.


__output file__

```
b 84600.000 447000.000 84699.998 447099.999
v 84693.637 447002.81 1.426
v 84693.613 447003.032 1.442
v 84693.468 447002.902 1.451
v 84693.564 447002.678 1.415
f 1 2 3
f 1 3 4
x 1
v 84692.714 447003.735 1.768
v 84692.743 447003.653 1.756
...
f 7 8 9
f 7 11 12
f 7 12 13
f 7 13 8
x 7
...
```

Vertices and triangles are 1-based; `x <id>` marks a vertex as finalised (it will no longer be referenced).


## 3. **ttsimpl**

  - optional stage between ttfin and ttdt
  - reads the .spa stream (stdin) and simplifies the points of each cell with a vertical epsilon
  - writes a reduced .spa stream to stdout (the other records are passed through unchanged)

Usage: `ttsimpl <vepsilon>`

```bash
./target/release/ttfin ./data/crop.laz 10 | ./target/release/ttsimpl 0.5 | ./target/release/ttdt > ./data/crop_10.sma
```

`<vepsilon>` is the maximum allowed vertical error (in the z units); points whose
error against the current triangulation is below it are discarded.


## 4. **ttstats**

  - shows some basic statistics to help you find the optimal cell size
  - reads a LAZ or XYZ file

```bash
./target/release/ttstats ./data/crop.laz
count:    212550
x-extent: 99.998
y-extent: 99.999
| cells | resolution | avg/cell | max/cell |
   8x8       13m          3321       4875
  16x16       7m           966       1676
  32x32       4m           363        661
  64x64       2m            92        224
 128x128      1m            23         87
 256x256      1m            23         87
 512x512      1m            23         87  
```


## 5. **ttval**

  - validates a .sma stream (stdin), eg to check that no vertex is used after it was finalised
  - prints the number of vertices and triangles

```bash
./target/release/ttfin ./data/crop.laz 10 | ./target/release/ttdt | ./target/release/ttval -vv
```
