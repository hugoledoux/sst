# sst

streaming startin



```bash
export RUST_LOG=info
./target/release/sstfin ../../data/ahn3/crop.laz 20 | ./target/release/sstdt | ./target/release/sstobj > ~/temp/hugo.obj
```


```bash
time RUST_LOG=info ./target/release/sstfin ~/data/ahn3/crop.txt 20 | RUST_LOG=info ./target/release/sstdt > ~/temp/z.txt
```

```bash
RUST_LOG=info ./target/debug/sstfin myinput.txt 10 | RUST_LOG=info ./target/debug/sstdel
```

