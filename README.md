# sst

__s__treaming __st__artin


## To compile them all

```bash
cargo build --release
```

`--release` makes the whole pipeline a whoooole lot faster, but you can't debug


## Four different binaries

### 1. __sstfin__

  - equivalent to spfinalise
  - finalises a set of points in LAZ or XYZ
  - can read several LAZ files at the same time
  - XYZ file is a CSV without header, space separated (see `/data/s400.xyz` for a simple example with 400 points)
  - reads the files on disk

To see the format of the output stream ("spa" is __s__treaming __p__oint __a__scii; binary to follow at some point)

```bash
export RUST_LOG=info
./target/release/sstfin ./data/crop.laz 10 > ./data/crop_10.spa
```

### 2. __sstdt__

  - equivalent to spdelaunay
  - takes as input (stdin) the .spa output of sstfin and creates a DT of the points
  - outputs a .sma (__s__treaming __m__esh __a__scii)

```bash
export RUST_LOG=info
./target/release/sstfin ./data/crop.laz 10 |  ./target/release/sstdt > ./data/crop_10.sma
```

### 3. __sstats__

  - shows some basics statistics to help you find the optimal cell size
  - reads a LAZ or XYZ file

```bash
./target/release/sstats ./data/crop.laz
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

### 4. __sstobj__

  - takes as input (stdin) the .sma output of sstdt and creates an OBJ file    

```bash
export RUST_LOG=info
./target/release/sstfin ./data/crop.laz 10 |  ./target/release/sstdt | ./target/release/sstobj > ~/temp/crop.obj
```

