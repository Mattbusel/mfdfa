//! Multifractal detrended fluctuation analysis (Kantelhardt et al. 2002).
//!
//! Given a stationary-increment signal `x` (for prices: log returns), MFDFA
//!
//! 1. builds the profile `Y(i) = sum_{k<=i} (x_k - mean(x))`;
//! 2. cuts `Y` into `N_s = floor(N / s)` non-overlapping segments of length
//!    `s`, once from the start and once from the end (`2 N_s` segments);
//! 3. removes a least-squares polynomial of degree `order` from each segment
//!    and keeps the residual variance `F^2(v, s)`;
//! 4. averages `F_q(s) = { mean_v [F^2(v, s)]^(q/2) }^(1/q)`, with the
//!    logarithmic average `F_0(s) = exp{ mean_v ln F^2(v, s) / 2 }` at `q = 0`;
//! 5. fits `h(q)` as the slope of `ln F_q(s)` against `ln s`.
//!
//! From `h(q)` the mass exponent is `tau(q) = q h(q) - 1` and the singularity
//! spectrum follows from the Legendre transform `alpha = d tau / d q`,
//! `f(alpha) = q alpha - tau(q)`. The derivative uses the same second-order
//! finite differences as `numpy.gradient`, so `alpha` and `f_alpha` have one
//! entry per `q`.
//!
//! `F_q(s)` matches the `MFDFA` Python package (Rydin Gorjao et al.) to
//! floating point precision; see `tests/reference.rs`.

use crate::error::{check_finite, ols_slope, Error, Result};

/// Output of [`mfdfa`].
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MfdfaSpectrum {
    /// The q orders, as passed in.
    pub q: Vec<f64>,
    /// The segment sizes, as passed in.
    pub scales: Vec<usize>,
    /// Fluctuation function: `fq[i][j] = F_{q_i}(s_j)`.
    pub fq: Vec<Vec<f64>>,
    /// Generalised Hurst exponents `h(q)`.
    pub h_q: Vec<f64>,
    /// Mass exponents `tau(q) = q h(q) - 1`.
    pub tau_q: Vec<f64>,
    /// Singularity strengths `alpha(q) = d tau / d q`.
    pub alpha: Vec<f64>,
    /// Singularity spectrum `f(alpha) = q alpha - tau(q)`.
    pub f_alpha: Vec<f64>,
}

impl MfdfaSpectrum {
    /// Spectrum width `max(alpha) - min(alpha)`. Near 0 for a monofractal
    /// signal; wider means more heterogeneous scaling.
    pub fn width(&self) -> f64 {
        let lo = self.alpha.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = self.alpha.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        hi - lo
    }

    /// `h(q)` at the q closest to 2, which is the ordinary DFA Hurst exponent.
    pub fn hurst(&self) -> f64 {
        let mut best = 0;
        for (i, q) in self.q.iter().enumerate() {
            if (q - 2.0).abs() < (self.q[best] - 2.0).abs() {
                best = i;
            }
        }
        self.h_q[best]
    }

    /// `h(q_min) - h(q_max)`: the spread of the generalised Hurst exponents.
    pub fn delta_h(&self) -> f64 {
        self.h_q[0] - self.h_q[self.h_q.len() - 1]
    }

    /// `alpha` at the peak of `f(alpha)` (the most common singularity).
    pub fn alpha_peak(&self) -> f64 {
        let mut best = 0;
        for (i, f) in self.f_alpha.iter().enumerate() {
            if *f > self.f_alpha[best] {
                best = i;
            }
        }
        self.alpha[best]
    }
}

/// Integrated profile `Y(i) = cumsum(x - mean(x))`.
pub fn profile(x: &[f64]) -> Vec<f64> {
    if x.is_empty() {
        return Vec::new();
    }
    let mean = x.iter().sum::<f64>() / x.len() as f64;
    let mut acc = 0.0;
    x.iter()
        .map(|v| {
            acc += v - mean;
            acc
        })
        .collect()
}

/// The default q grid: -5 to 5 in steps of 0.5 (21 values).
pub fn default_q() -> Vec<f64> {
    (0..=20).map(|i| -5.0 + 0.5 * i as f64).collect()
}

/// Up to 20 log-spaced scales from `max(10, 2 (order + 1))` to `n / 4`,
/// rounded and de-duplicated. May return fewer than two scales for very short
/// series, in which case [`mfdfa`] reports an error.
pub fn default_scales(n: usize, order: usize) -> Vec<usize> {
    let lo = 10usize.max(2 * (order + 1));
    let hi = n / 4;
    if hi < lo {
        return Vec::new();
    }
    let (llo, lhi) = ((lo as f64).ln(), (hi as f64).ln());
    let points = 20;
    let mut out: Vec<usize> = (0..points)
        .map(|i| {
            (llo + (lhi - llo) * i as f64 / (points - 1) as f64)
                .exp()
                .round() as usize
        })
        .collect();
    out.dedup();
    out
}

/// Orthonormal basis (columns) of polynomials of degree `0..=order` sampled on
/// `s` points, built with twice-repeated modified Gram-Schmidt on a scaled
/// abscissa for numerical stability.
fn poly_basis(s: usize, order: usize) -> Vec<Vec<f64>> {
    let t: Vec<f64> = if s > 1 {
        (0..s)
            .map(|i| 2.0 * i as f64 / (s - 1) as f64 - 1.0)
            .collect()
    } else {
        vec![0.0]
    };
    let mut basis: Vec<Vec<f64>> = Vec::with_capacity(order + 1);
    for d in 0..=order {
        let mut v: Vec<f64> = t.iter().map(|ti| ti.powi(d as i32)).collect();
        for _ in 0..2 {
            for b in &basis {
                let c: f64 = v.iter().zip(b).map(|(a, b)| a * b).sum();
                for (vi, bi) in v.iter_mut().zip(b) {
                    *vi -= c * bi;
                }
            }
        }
        let norm = v.iter().map(|a| a * a).sum::<f64>().sqrt();
        for vi in v.iter_mut() {
            *vi /= norm;
        }
        basis.push(v);
    }
    basis
}

/// Residual variance of `seg` after removing its least-squares polynomial.
fn detrended_variance(seg: &[f64], basis: &[Vec<f64>], resid: &mut [f64]) -> f64 {
    resid.copy_from_slice(seg);
    for b in basis {
        let c: f64 = resid.iter().zip(b).map(|(a, b)| a * b).sum();
        for (r, bi) in resid.iter_mut().zip(b) {
            *r -= c * bi;
        }
    }
    resid.iter().map(|r| r * r).sum::<f64>() / seg.len() as f64
}

/// Segment variances `F^2(v, s)` for one scale (forward then backward).
fn segment_variances(prof: &[f64], s: usize, order: usize) -> Vec<f64> {
    let n = prof.len();
    let ns = n / s;
    let basis = poly_basis(s, order);
    let mut resid = vec![0.0; s];
    let mut out = Vec::with_capacity(2 * ns);
    for v in 0..ns {
        out.push(detrended_variance(
            &prof[v * s..(v + 1) * s],
            &basis,
            &mut resid,
        ));
    }
    let off = n % s;
    for v in 0..ns {
        let start = off + v * s;
        out.push(detrended_variance(
            &prof[start..start + s],
            &basis,
            &mut resid,
        ));
    }
    out
}

fn q_average(f2: &[f64], q: f64) -> f64 {
    let m = f2.len() as f64;
    if q == 0.0 {
        (0.5 * f2.iter().map(|v| v.ln()).sum::<f64>() / m).exp()
    } else {
        (f2.iter().map(|v| v.powf(q / 2.0)).sum::<f64>() / m).powf(1.0 / q)
    }
}

fn validate(x: &[f64], q: &[f64], scales: &[usize], order: usize) -> Result<()> {
    check_finite(x)?;
    if q.is_empty() {
        return Err(Error::InvalidParameter("q must not be empty"));
    }
    if q.iter().any(|v| !v.is_finite()) {
        return Err(Error::InvalidParameter("q values must be finite"));
    }
    if scales.len() < 2 {
        return Err(Error::InvalidParameter("need at least two scales"));
    }
    let mut sorted = scales.to_vec();
    sorted.sort_unstable();
    if sorted.windows(2).any(|w| w[0] == w[1]) {
        return Err(Error::InvalidParameter("scales must be distinct"));
    }
    let min_s = order + 2;
    let needed = 2 * sorted[sorted.len() - 1];
    if sorted[0] < min_s {
        return Err(Error::InvalidParameter(
            "every scale must exceed the detrending order by at least 2",
        ));
    }
    if x.len() < needed {
        return Err(Error::TooShort {
            needed,
            got: x.len(),
        });
    }
    Ok(())
}

/// Fluctuation function `F_q(s)` for every `q` (rows) and scale (columns).
pub fn fluctuation_function(
    x: &[f64],
    q: &[f64],
    scales: &[usize],
    order: usize,
) -> Result<Vec<Vec<f64>>> {
    validate(x, q, scales, order)?;
    let prof = profile(x);
    let mut fq = vec![vec![0.0; scales.len()]; q.len()];
    for (j, &s) in scales.iter().enumerate() {
        let f2 = segment_variances(&prof, s, order);
        for (i, &qi) in q.iter().enumerate() {
            let v = q_average(&f2, qi);
            if !(v.is_finite() && v > 0.0) {
                return Err(Error::Degenerate(
                    "zero or non-finite fluctuation (constant segment with q <= 0?)",
                ));
            }
            fq[i][j] = v;
        }
    }
    Ok(fq)
}

/// `numpy.gradient(f, x)` with first-order edges.
pub(crate) fn gradient(f: &[f64], x: &[f64]) -> Vec<f64> {
    let n = f.len();
    if n < 2 {
        return vec![f64::NAN; n];
    }
    let mut g = vec![0.0; n];
    g[0] = (f[1] - f[0]) / (x[1] - x[0]);
    g[n - 1] = (f[n - 1] - f[n - 2]) / (x[n - 1] - x[n - 2]);
    for i in 1..n - 1 {
        let hd = x[i] - x[i - 1];
        let hs = x[i + 1] - x[i];
        g[i] = (hd * hd * f[i + 1] + (hs * hs - hd * hd) * f[i] - hs * hs * f[i - 1])
            / (hs * hd * (hd + hs));
    }
    g
}

/// Multifractal detrended fluctuation analysis.
///
/// * `x`: the signal (increments; for prices pass log returns).
/// * `q`: q orders, strictly increasing, at least two (for example
///   [`default_q`]). `q = 0` is handled with the logarithmic average.
/// * `scales`: distinct segment sizes, each at least `order + 2` and at most
///   `x.len() / 2` (for example [`default_scales`]).
/// * `order`: degree of the detrending polynomial (1 = MFDFA-1).
///
/// ```
/// use mfdfa::{mfdfa, default_q, default_scales};
/// // Deterministic white noise from a linear congruential generator.
/// let mut s = 1u64;
/// let x: Vec<f64> = (0..2048)
///     .map(|_| {
///         s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
///         (s >> 11) as f64 / (1u64 << 53) as f64 - 0.5
///     })
///     .collect();
/// let spec = mfdfa(&x, &default_q(), &default_scales(x.len(), 1), 1).unwrap();
/// assert_eq!(spec.h_q.len(), 21);
/// assert!((spec.hurst() - 0.5).abs() < 0.1);
/// ```
pub fn mfdfa(x: &[f64], q: &[f64], scales: &[usize], order: usize) -> Result<MfdfaSpectrum> {
    if q.len() < 2 {
        return Err(Error::InvalidParameter(
            "need at least two q values for the spectrum",
        ));
    }
    if q.windows(2).any(|w| w[1] <= w[0]) {
        return Err(Error::InvalidParameter("q must be strictly increasing"));
    }
    let fq = fluctuation_function(x, q, scales, order)?;
    let log_s: Vec<f64> = scales.iter().map(|&s| (s as f64).ln()).collect();
    let h_q: Vec<f64> = fq
        .iter()
        .map(|row| {
            let ly: Vec<f64> = row.iter().map(|v| v.ln()).collect();
            ols_slope(&log_s, &ly)
        })
        .collect();
    let tau_q: Vec<f64> = q.iter().zip(&h_q).map(|(q, h)| q * h - 1.0).collect();
    let alpha = gradient(&tau_q, q);
    let f_alpha: Vec<f64> = q
        .iter()
        .zip(&alpha)
        .zip(&tau_q)
        .map(|((q, a), t)| q * a - t)
        .collect();
    Ok(MfdfaSpectrum {
        q: q.to_vec(),
        scales: scales.to_vec(),
        fq,
        h_q,
        tau_q,
        alpha,
        f_alpha,
    })
}

/// Ordinary DFA Hurst exponent: the slope of `ln F_2(s)` on `ln s`.
///
/// `x` is the increment series (returns). An exponent near 0.5 means no
/// long memory, above 0.5 persistence, below 0.5 anti-persistence.
pub fn dfa_hurst(x: &[f64], scales: &[usize], order: usize) -> Result<f64> {
    let fq = fluctuation_function(x, &[2.0], scales, order)?;
    let log_s: Vec<f64> = scales.iter().map(|&s| (s as f64).ln()).collect();
    let ly: Vec<f64> = fq[0].iter().map(|v| v.ln()).collect();
    Ok(ols_slope(&log_s, &ly))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg_noise(n: usize, seed: u64) -> Vec<f64> {
        let mut s = seed;
        (0..n)
            .map(|_| {
                s = s
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                (s >> 11) as f64 / (1u64 << 53) as f64 - 0.5
            })
            .collect()
    }

    #[test]
    fn profile_ends_at_zero() {
        let p = profile(&lcg_noise(500, 3));
        assert!(p[p.len() - 1].abs() < 1e-10);
    }

    #[test]
    fn basis_is_orthonormal() {
        let b = poly_basis(17, 3);
        for i in 0..4 {
            for j in 0..4 {
                let d: f64 = b[i].iter().zip(&b[j]).map(|(a, c)| a * c).sum();
                let want = if i == j { 1.0 } else { 0.0 };
                assert!((d - want).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn polynomial_is_removed_exactly() {
        let seg: Vec<f64> = (0..20)
            .map(|i| 3.0 - 0.5 * i as f64 + 0.01 * (i * i) as f64)
            .collect();
        let b = poly_basis(20, 2);
        let mut r = vec![0.0; 20];
        assert!(detrended_variance(&seg, &b, &mut r) < 1e-20);
    }

    #[test]
    fn gradient_matches_numpy_formula() {
        let x = [-2.0, -1.0, 0.0, 0.5, 2.0];
        let f: Vec<f64> = x.iter().map(|v: &f64| v * v).collect();
        let g = gradient(&f, &x);
        // Second-order interior differences are exact for a quadratic.
        assert!((g[1] - -2.0).abs() < 1e-12);
        assert!((g[2] - 0.0).abs() < 1e-12);
        assert!((g[3] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn white_noise_is_monofractal_near_half() {
        let x = lcg_noise(8192, 11);
        let spec = mfdfa(&x, &default_q(), &default_scales(x.len(), 1), 1).unwrap();
        assert!((spec.hurst() - 0.5).abs() < 0.06, "h(2) = {}", spec.hurst());
        assert!(spec.width() < 0.35, "width = {}", spec.width());
    }

    #[test]
    fn rejects_bad_parameters() {
        let x = lcg_noise(100, 1);
        assert!(matches!(
            mfdfa(&x, &[2.0], &[10, 20], 1),
            Err(Error::InvalidParameter(_))
        ));
        assert!(matches!(
            mfdfa(&x, &[2.0, 1.0], &[10, 20], 1),
            Err(Error::InvalidParameter(_))
        ));
        assert!(matches!(
            mfdfa(&x, &[1.0, 2.0], &[10, 10], 1),
            Err(Error::InvalidParameter(_))
        ));
        assert!(matches!(
            mfdfa(&x, &[1.0, 2.0], &[2, 10], 1),
            Err(Error::InvalidParameter(_))
        ));
        assert!(matches!(
            mfdfa(&x, &[1.0, 2.0], &[10, 60], 1),
            Err(Error::TooShort { .. })
        ));
        let mut bad = x.clone();
        bad[5] = f64::NAN;
        assert_eq!(
            mfdfa(&bad, &[1.0, 2.0], &[10, 20], 1),
            Err(Error::NonFiniteInput)
        );
    }

    #[test]
    fn constant_series_is_degenerate_for_negative_q() {
        let x = vec![1.0; 200];
        assert!(matches!(
            mfdfa(&x, &[-2.0, 2.0], &[10, 20], 1),
            Err(Error::Degenerate(_))
        ));
    }
}
