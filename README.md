# sst

**s**treaming **st**artin


## To compile them all

```bash
cargo build --release
```

`--release` makes the whole pipeline a whoooole lot faster, but you can't debug


## Four different binaries

### 1. **sstfin**

  - equivalent to spfinalise
  - finalises a set of points in LAZ or XYZ
  - can read several LAZ files at the same time
  - XYZ file is a CSV without header, space separated (see `/data/s400.xyz` for a simple example with 400 points)
  - reads the files on disk

To see the format of the output stream ("spa" is **s**treaming **p**oint **a**scii; binary to follow at some point)

```bash
export RUST_LOG=info
./target/release/sstfin ./data/crop.laz 10 > ./data/crop_10.spa
```

Multiple LAZ files if you provide them in a `myfiles.files`, see the example in `./data/crop_1_2.files`.
This reads `crop_1.laz` and then `crop_2.laz` (in that order); I split `crop.laz` to obtain those 2 files.

```bash
export RUST_LOG=info
./target/release/sstfin ./data/crop_1_2.files 10 > /dev/null/
```

### 2. **sstdt**

  - equivalent to spdelaunay
  - takes as input (stdin) the .spa output of sstfin and creates a DT of the points
  - outputs a .sma (**s**treaming **m**esh **a**scii)

```bash
export RUST_LOG=info
./target/release/sstfin ./data/crop.laz 10 |  ./target/release/sstdt > ./data/crop_10.sma
```

### 3. **sstats**

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

### 4. **sstobj**

  - takes as input (stdin) the .sma output of sstdt and creates an OBJ file    

```bash
export RUST_LOG=info
./target/release/sstfin ./data/crop.laz 10 |  ./target/release/sstdt | ./target/release/sstobj > ~/temp/crop.obj
```

