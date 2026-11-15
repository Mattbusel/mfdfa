//! # mfdfa
//!
//! Multifractal detrended fluctuation analysis, Higuchi fractal dimension and
//! wavelet Hurst estimation for time series, with no required dependencies.
//!
//! | Function | What it measures | Input |
//! |---|---|---|
//! | [`mfdfa`] | generalised Hurst `h(q)`, `tau(q)`, singularity spectrum `f(alpha)` | increments (returns) |
//! | [`dfa_hurst`] | ordinary DFA Hurst exponent (`h(2)`) | increments |
//! | [`higuchi_fd`] / [`higuchi`] | Higuchi fractal dimension of a curve | path (prices) |
//! | [`wavelet_hurst`] | Abry-Veitch wavelet Hurst exponent (db4) | increments |
//! | [`hurst::hurst_rs`] | rescaled-range Hurst exponent | prices |
//! | [`RollingHurst`], [`RollingHiguchi`] | streaming versions (`ta` feature adds `ta::Next`) | bar by bar |
//!
//! ```
//! use mfdfa::{mfdfa, default_q, default_scales, higuchi_fd};
//!
//! let mut s = 42u64;
//! let returns: Vec<f64> = (0..4096)
//!     .map(|_| {
//!         s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
//!         ((s >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * 0.02
//!     })
//!     .collect();
//!
//! let spec = mfdfa(&returns, &default_q(), &default_scales(returns.len(), 1), 1).unwrap();
//! println!("h(2) = {:.3}, width = {:.3}", spec.hurst(), spec.width());
//!
//! // Higuchi works on the path, not the increments.
//! let mut p = 0.0;
//! let path: Vec<f64> = returns.iter().map(|r| { p += r; p }).collect();
//! let d = higuchi_fd(&path, 10);
//! assert!((d - 1.5).abs() < 0.1); // Brownian path
//! ```
//!
//! ## Features
//!
//! * `ta`: `ta::Next<f64>`, `ta::Next<&T: Close>`, `ta::Reset`, `ta::Period`
//!   for the rolling estimators.
//! * `serde`: `Serialize` / `Deserialize` on result types.
//! * `synth`: exact fractional Gaussian noise (`synth::fgn`, `synth::fbm`).
//! * `cli`: the `mfdfa` command line tool.
//!
//! ## Provenance
//!
//! Carved out of the `fractal-engine` crate and `lib/math/fractal.py` in
//! [srfm-lab](https://github.com/Mattbusel/srfm-lab) (MIT). Numerical results
//! are checked against the `MFDFA`, `antropy` and `PyWavelets` Python
//! packages in `tests/reference.rs`.

#![cfg_attr(docsrs, feature(doc_cfg))]
#![warn(missing_docs)]

mod dfa;
mod error;
mod higuchi;
mod streaming;

pub mod detrending;
pub mod fractal_dim_ext;
pub mod fractal_dimension;
pub mod hurst;
pub mod multifractal;
#[cfg(feature = "synth")]
#[cfg_attr(docsrs, doc(cfg(feature = "synth")))]
pub mod synth;
pub mod wavelet;

pub use dfa::{
    default_q, default_scales, dfa_hurst, fluctuation_function, mfdfa, profile, MfdfaSpectrum,
};
pub use error::{Error, Result};
pub use higuchi::{higuchi, higuchi_fd, Higuchi};
pub use streaming::{Input, RollingHiguchi, RollingHurst};
pub use wavelet::{wavelet_hurst, wavelet_hurst_detail, WaveletHurst};
