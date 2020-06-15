# sst

streaming startin

```bash
RUST_LOG=info ./target/debug/sstfin ./data/square400.txt 15 | RUST_LOG=info ./target/debug/sstdt | RUST_LOG=info ./target/debug/sstobj > ~/temp/o.obj
```


```bash
time RUST_LOG=info ./target/release/sstfin ~/data/ahn3/crop.txt 20 | RUST_LOG=info ./target/release/sstdel > ~/temp/z.txt
```

```bash
RUST_LOG=info ./target/debug/sstfin myinput.txt 10 | RUST_LOG=info ./target/debug/sstdel
```

