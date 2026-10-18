//! Binary variants of the two stream formats (`.spa` and `.sma`).
//!
//! A binary stream is a sequence of *fixed-size little-endian records*:
//! one tag byte followed by a fixed payload, preceded by a magic header
//! (`TTSPA1\n` / `TTSMA1\n`). The records are the exact equivalents of the
//! ASCII lines (see the `to_line()` methods), so a reader can auto-detect
//! the format by peeking the first bytes and a tool that only understands
//! the ASCII lines can still consume a binary stream through
//! `record.to_line()` (this is what the viewers do).
//!
//! Writers pick the encoding with their `--ascii` flag (binary is the
//! default); readers auto-detect and need no flag.

use std::io::{self, BufRead, Write};

pub const SPA_MAGIC: &[u8] = b"TTSPA1\n";
pub const SMA_MAGIC: &[u8] = b"TTSMA1\n";

//-- .spa records -------------------------------------------------------------
//-- tag + payload (little-endian):
//--   'n' u64          'c' u32          's' u64
//--   'b' 4 x f64      'v' 3 x f64      'w' 3 x f64      'x' 2 x u32

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SpaRecord {
    N(u64),
    C(u32),
    S(u64),
    B([f64; 4]),
    V(f64, f64, f64),
    W(f64, f64, f64),
    X(u32, u32),
}

impl SpaRecord {
    /// The same record as an ASCII line, as written by `ttfin --ascii`.
    pub fn to_line(&self) -> String {
        match self {
            SpaRecord::N(n) => format!("n {}\n", n),
            SpaRecord::C(c) => format!("c {}\n", c),
            SpaRecord::S(s) => format!("s {}\n", s),
            SpaRecord::B(b) => format!("b {:.3} {:.3} {:.3} {:.3}\n", b[0], b[1], b[2], b[3]),
            SpaRecord::V(x, y, z) => format!("v {} {} {}\n", x, y, z),
            SpaRecord::W(x, y, z) => format!("w {} {} {}\n", x, y, z),
            SpaRecord::X(a, b) => format!("x {} {}\n", a, b),
        }
    }
}

//-- .sma records ------------------------------------------------------------
//-- tag + payload (little-endian):
//--   'b' 4 x f64      'v' 3 x f64      'f' 3 x u32      'x' u32

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SmaRecord {
    B([f64; 4]),
    V(f64, f64, f64),
    F(u32, u32, u32),
    X(u32),
}

impl SmaRecord {
    /// The same record as an ASCII line, as written by `ttdt --ascii`.
    pub fn to_line(&self) -> String {
        match self {
            SmaRecord::B(b) => format!("b {:.3} {:.3} {:.3} {:.3}\n", b[0], b[1], b[2], b[3]),
            SmaRecord::V(x, y, z) => format!("v {} {} {}\n", x, y, z),
            SmaRecord::F(a, b, c) => format!("f {} {} {}\n", a, b, c),
            SmaRecord::X(a) => format!("x {}\n", a),
        }
    }
}

//-- writing ------------------------------------------------------------------

const MAX_REC: usize = 33; //-- largest record: tag + 4 x f64

pub fn spa_write_magic<W: Write>(w: &mut W) -> io::Result<()> {
    w.write_all(SPA_MAGIC)
}

pub fn sma_write_magic<W: Write>(w: &mut W) -> io::Result<()> {
    w.write_all(SMA_MAGIC)
}

pub fn spa_write<W: Write>(w: &mut W, r: &SpaRecord) -> io::Result<()> {
    let mut buf = [0u8; MAX_REC];
    let n = match r {
        SpaRecord::N(v) => {
            buf[0] = b'n';
            buf[1..9].copy_from_slice(&v.to_le_bytes());
            9
        }
        SpaRecord::C(v) => {
            buf[0] = b'c';
            buf[1..5].copy_from_slice(&v.to_le_bytes());
            5
        }
        SpaRecord::S(v) => {
            buf[0] = b's';
            buf[1..9].copy_from_slice(&v.to_le_bytes());
            9
        }
        SpaRecord::B(b) => {
            buf[0] = b'b';
            for i in 0..4 {
                buf[1 + 8 * i..9 + 8 * i].copy_from_slice(&b[i].to_le_bytes());
            }
            33
        }
        SpaRecord::V(x, y, z) => {
            buf[0] = b'v';
            buf[1..9].copy_from_slice(&x.to_le_bytes());
            buf[9..17].copy_from_slice(&y.to_le_bytes());
            buf[17..25].copy_from_slice(&z.to_le_bytes());
            25
        }
        SpaRecord::W(x, y, z) => {
            buf[0] = b'w';
            buf[1..9].copy_from_slice(&x.to_le_bytes());
            buf[9..17].copy_from_slice(&y.to_le_bytes());
            buf[17..25].copy_from_slice(&z.to_le_bytes());
            25
        }
        SpaRecord::X(a, b) => {
            buf[0] = b'x';
            buf[1..5].copy_from_slice(&a.to_le_bytes());
            buf[5..9].copy_from_slice(&b.to_le_bytes());
            9
        }
    };
    w.write_all(&buf[..n])
}

pub fn sma_write<W: Write>(w: &mut W, r: &SmaRecord) -> io::Result<()> {
    let mut buf = [0u8; MAX_REC];
    let n = match r {
        SmaRecord::B(b) => {
            buf[0] = b'b';
            for i in 0..4 {
                buf[1 + 8 * i..9 + 8 * i].copy_from_slice(&b[i].to_le_bytes());
            }
            33
        }
        SmaRecord::V(x, y, z) => {
            buf[0] = b'v';
            buf[1..9].copy_from_slice(&x.to_le_bytes());
            buf[9..17].copy_from_slice(&y.to_le_bytes());
            buf[17..25].copy_from_slice(&z.to_le_bytes());
            25
        }
        SmaRecord::F(a, b, c) => {
            buf[0] = b'f';
            buf[1..5].copy_from_slice(&a.to_le_bytes());
            buf[5..9].copy_from_slice(&b.to_le_bytes());
            buf[9..13].copy_from_slice(&c.to_le_bytes());
            13
        }
        SmaRecord::X(a) => {
            buf[0] = b'x';
            buf[1..5].copy_from_slice(&a.to_le_bytes());
            5
        }
    };
    w.write_all(&buf[..n])
}

//-- reading ------------------------------------------------------------------

/// Tell binary from ASCII without consuming the buffer: peek the first
/// bytes of a buffered source and compare with the magic of the format the
/// caller expects (`SPA_MAGIC` or `SMA_MAGIC`).
pub fn peek_is_binary<R: BufRead>(r: &mut R, magic: &[u8]) -> io::Result<bool> {
    let buf = r.fill_buf()?;
    Ok(buf.len() >= magic.len() && &buf[..magic.len()] == magic)
}

/// Consume the magic (only after `peek_is_binary` said yes).
pub fn skip_magic<R: BufRead>(r: &mut R, magic: &[u8]) -> io::Result<()> {
    let mut buf = vec![0u8; magic.len()];
    r.read_exact(&mut buf)?;
    Ok(())
}

fn read_byte<R: BufRead>(r: &mut R) -> io::Result<Option<u8>> {
    let mut b = [0u8; 1];
    Ok(match r.read(&mut b)? {
        0 => None,
        _ => Some(b[0]),
    })
}

fn read_u32<R: BufRead>(r: &mut R) -> io::Result<u32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}

fn read_u64<R: BufRead>(r: &mut R) -> io::Result<u64> {
    let mut b = [0u8; 8];
    r.read_exact(&mut b)?;
    Ok(u64::from_le_bytes(b))
}

fn read_f64<R: BufRead>(r: &mut R) -> io::Result<f64> {
    let mut b = [0u8; 8];
    r.read_exact(&mut b)?;
    Ok(f64::from_le_bytes(b))
}

fn read_3f64<R: BufRead>(r: &mut R) -> io::Result<(f64, f64, f64)> {
    Ok((read_f64(r)?, read_f64(r)?, read_f64(r)?))
}

fn read_4f64<R: BufRead>(r: &mut R) -> io::Result<[f64; 4]> {
    Ok([
        read_f64(r)?,
        read_f64(r)?,
        read_f64(r)?,
        read_f64(r)?,
    ])
}

/// Read one binary `.spa` record; `None` on a clean EOF at a record boundary.
/// A truncated record surfaces as an `UnexpectedEof` error.
pub fn spa_read<R: BufRead>(r: &mut R) -> io::Result<Option<SpaRecord>> {
    let tag = match read_byte(r)? {
        None => return Ok(None),
        Some(b) => b,
    };
    let rec = match tag {
        b'n' => SpaRecord::N(read_u64(r)?),
        b'c' => SpaRecord::C(read_u32(r)?),
        b's' => SpaRecord::S(read_u64(r)?),
        b'b' => SpaRecord::B(read_4f64(r)?),
        b'v' => {
            let (x, y, z) = read_3f64(r)?;
            SpaRecord::V(x, y, z)
        }
        b'w' => {
            let (x, y, z) = read_3f64(r)?;
            SpaRecord::W(x, y, z)
        }
        b'x' => {
            let a = read_u32(r)?;
            let b = read_u32(r)?;
            SpaRecord::X(a, b)
        }
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("unknown .spa record tag '{}'", tag as char),
            ))
        }
    };
    Ok(Some(rec))
}

/// Read one binary `.sma` record; `None` on a clean EOF at a record boundary.
pub fn sma_read<R: BufRead>(r: &mut R) -> io::Result<Option<SmaRecord>> {
    let tag = match read_byte(r)? {
        None => return Ok(None),
        Some(b) => b,
    };
    let rec = match tag {
        b'b' => SmaRecord::B(read_4f64(r)?),
        b'v' => {
            let (x, y, z) = read_3f64(r)?;
            SmaRecord::V(x, y, z)
        }
        b'f' => {
            let a = read_u32(r)?;
            let b = read_u32(r)?;
            let c = read_u32(r)?;
            SmaRecord::F(a, b, c)
        }
        b'x' => SmaRecord::X(read_u32(r)?),
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("unknown .sma record tag '{}'", tag as char),
            ))
        }
    };
    Ok(Some(rec))
}
