//! Higuchi fractal dimension (Higuchi 1988, Physica D 31, 277-283).
//!
//! For each delay `k = 1..=k_max` and offset `m = 1..=k` the curve length is
//!
//! ```text
//! L_m(k) = (N - 1) / (k^2 floor((N - m) / k)) * sum_{j=1}^{floor((N-m)/k)} |x(m + j k) - x(m + (j-1) k)|
//! ```
//!
//! `L(k)` is the mean of `L_m(k)` over `m`, and the fractal dimension is the
//! slope of `ln L(k)` against `ln(1 / k)`. The input is the curve itself
//! (for example a price or log-price path, not its returns). A smooth line
//! gives D close to 1, Brownian motion about 1.5 and white noise about 2;
//! for fractional Brownian motion D is about `2 - H`.
//!
//! Matches `antropy.higuchi_fd` to floating point precision (see
//! `tests/reference.rs`).

use crate::error::{check_finite, ols_slope, Error, Result};

/// Detailed output of [`higuchi`].
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Higuchi {
    /// Estimated fractal dimension (not clamped).
    pub fd: f64,
    /// Delays used, `1..=k_max`.
    pub k: Vec<usize>,
    /// Mean curve length `L(k)` for each delay.
    pub curve_length: Vec<f64>,
}

/// Mean normalised curve length `L(k)`.
pub(crate) fn curve_length(x: &[f64], k: usize) -> f64 {
    let n = x.len();
    let mut total = 0.0;
    for m in 0..k {
        let steps = (n - m - 1) / k;
        let mut l = 0.0;
        for j in 1..=steps {
            l += (x[m + j * k] - x[m + (j - 1) * k]).abs();
        }
        l *= (n - 1) as f64 / (k * k * steps) as f64;
        total += l;
    }
    total / k as f64
}

/// Higuchi fractal dimension with the curve lengths behind it.
///
/// Needs `k_max >= 2` and `x.len() >= 2 * k_max` so that every offset has at
/// least one increment.
pub fn higuchi(x: &[f64], k_max: usize) -> Result<Higuchi> {
    check_finite(x)?;
    if k_max < 2 {
        return Err(Error::InvalidParameter("k_max must be at least 2"));
    }
    let needed = 2 * k_max;
    if x.len() < needed {
        return Err(Error::TooShort {
            needed,
            got: x.len(),
        });
    }
    let k: Vec<usize> = (1..=k_max).collect();
    let curve: Vec<f64> = k.iter().map(|&k| curve_length(x, k)).collect();
    if curve.iter().any(|l| l.is_nan() || *l <= 0.0) {
        return Err(Error::Degenerate("zero curve length (constant series)"));
    }
    let lx: Vec<f64> = k.iter().map(|&k| -(k as f64).ln()).collect();
    let ly: Vec<f64> = curve.iter().map(|l| l.ln()).collect();
    Ok(Higuchi {
        fd: ols_slope(&lx, &ly),
        k,
        curve_length: curve,
    })
}

/// Higuchi fractal dimension of the curve `x`, or NaN when it cannot be
/// computed (series shorter than `2 * k_max`, `k_max < 2`, constant or
/// non-finite input). Use [`higuchi`] for the error reason.
///
/// ```
/// let line: Vec<f64> = (0..200).map(|i| i as f64).collect();
/// assert!((mfdfa::higuchi_fd(&line, 10) - 1.0).abs() < 1e-9);
/// ```
pub fn higuchi_fd(x: &[f64], k_max: usize) -> f64 {
    higuchi(x, k_max).map(|h| h.fd).unwrap_or(f64::NAN)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn straight_line_is_one() {
        let x: Vec<f64> = (0..100).map(|i| 2.0 * i as f64 + 1.0).collect();
        assert!((higuchi_fd(&x, 8) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn white_noise_is_near_two() {
        let mut s = 7u64;
        let x: Vec<f64> = (0..4096)
            .map(|_| {
                s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
                (s >> 11) as f64 / (1u64 << 53) as f64
            })
            .collect();
        let fd = higuchi_fd(&x, 10);
        assert!((fd - 2.0).abs() < 0.05, "fd = {fd}");
    }

    #[test]
    fn errors() {
        assert_eq!(
            higuchi(&[1.0; 10], 1),
            Err(Error::InvalidParameter("k_max must be at least 2"))
        );
        assert_eq!(
            higuchi(&[1.0; 10], 8),
            Err(Error::TooShort {
                needed: 16,
                got: 10
            })
        );
        assert!(matches!(higuchi(&[1.0; 40], 8), Err(Error::Degenerate(_))));
        assert!(higuchi_fd(&[1.0; 40], 8).is_nan());
    }
}
