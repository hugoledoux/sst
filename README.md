# sst

**s**treaming **st**artin


## To compile them all

```bash
cargo build --release
```

`--release` makes the whole pipeline a whoooole lot faster, but you can't debug


## Three different binaries

### 1. **sstfin**

  - conceptually equivalent to Isenburg's *spfinalise*
  - finalises a set of points in LAZ or XYZ
  - can read several LAZ files at the same time
  - XYZ file is a CSV without header, space separated (see `/data/s400.xyz` for a simple example with 400 points)
  - reads the files on disk

To see the format of the output stream ("spa" is **s**treaming **p**oint **a**scii; binary to follow at some point)

```bash
./target/release/sstfin -vv ./data/crop.laz 10 > ./data/crop_10.spa
```

Multiple LAZ files if you provide them in a `myfiles.files`, see the example in `./data/crop_1_2.files`.
This reads `crop_1.laz` and then `crop_2.laz` (in that order); I split `crop.laz` to obtain those 2 files.

```bash
./target/release/sstfin -vv ./data/crop_1_2.files 10 > /dev/null
```

(remove `-vv` to remove the INFO and WARN)

__.spa file__

```
# sstfin
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

### 2. **sstdt**

  - conceptually equivalent to Isenburg's *spdelaunay*;
  - takes as input (stdin) the .spa output of sstfin and creates a DT of the points;
  - outputs a variation of an .obj file where a few vertex finalisers are added. It's simpler than Isenburg's .sma (**s**treaming **m**esh **a**scii) and can be opened as an .obj file by different software (MeshLab just ignores the extra lines, so you can open the file with it);

```bash
./target/release/sstfin -vv ./data/crop.laz 10 | ./target/release/sstdt -vv > ./data/crop_10.sma
```

The output file `crop_10.sma` is a valid OBJ file: change its extension and you can open it normally in [MeshLab](https://www.meshlab.net/) for instance.


__output file__

```
# sstdt
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

