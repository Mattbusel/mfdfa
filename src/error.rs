//! Error type shared by every estimator in the crate.

use core::fmt;

/// Why an estimate could not be produced.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// The input series is shorter than the estimator needs.
    TooShort {
        /// Minimum number of observations required.
        needed: usize,
        /// Number of observations supplied.
        got: usize,
    },
    /// A parameter is out of range (the message names it).
    InvalidParameter(&'static str),
    /// The input contains NaN or infinity.
    NonFiniteInput,
    /// An intermediate quantity was zero or non-finite, for example a
    /// constant series (zero fluctuation) combined with a negative `q`.
    Degenerate(&'static str),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::TooShort { needed, got } => {
                write!(
                    f,
                    "series too short: need at least {needed} observations, got {got}"
                )
            }
            Error::InvalidParameter(msg) => write!(f, "invalid parameter: {msg}"),
            Error::NonFiniteInput => write!(f, "input contains NaN or infinite values"),
            Error::Degenerate(msg) => write!(f, "degenerate input: {msg}"),
        }
    }
}

impl std::error::Error for Error {}

/// Crate result alias.
pub type Result<T> = core::result::Result<T, Error>;

pub(crate) fn check_finite(x: &[f64]) -> Result<()> {
    if x.iter().all(|v| v.is_finite()) {
        Ok(())
    } else {
        Err(Error::NonFiniteInput)
    }
}

/// Ordinary least squares slope of `y` on `x`.
pub(crate) fn ols_slope(x: &[f64], y: &[f64]) -> f64 {
    debug_assert_eq!(x.len(), y.len());
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let mut num = 0.0;
    let mut den = 0.0;
    for (a, b) in x.iter().zip(y) {
        num += (a - mx) * (b - my);
        den += (a - mx) * (a - mx);
    }
    num / den
}

/// Weighted least squares slope of `y` on `x` with weights `w`.
pub(crate) fn wls_slope(x: &[f64], y: &[f64], w: &[f64]) -> f64 {
    let sw: f64 = w.iter().sum();
    let mx = x.iter().zip(w).map(|(a, b)| a * b).sum::<f64>() / sw;
    let my = y.iter().zip(w).map(|(a, b)| a * b).sum::<f64>() / sw;
    let mut num = 0.0;
    let mut den = 0.0;
    for ((a, b), c) in x.iter().zip(y).zip(w) {
        num += c * (a - mx) * (b - my);
        den += c * (a - mx) * (a - mx);
    }
    num / den
}

/// Log returns of a positive price series; non-positive prices yield 0.
pub(crate) fn log_returns(prices: &[f64]) -> Vec<f64> {
    prices
        .windows(2)
        .map(|w| {
            if w[0] > 0.0 && w[1] > 0.0 {
                (w[1] / w[0]).ln()
            } else {
                0.0
            }
        })
        .collect()
}
