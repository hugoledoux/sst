//! Shared library code for the tintin tools.
//!
//! Currently the binary stream codecs (`stream`) shared by all the tools,
//! and the 2D viewer plumbing (`view`) shared by the standalone `ttview`
//! binary and the `ttdt --gui` embedded view.

pub mod stream;
pub mod view;
