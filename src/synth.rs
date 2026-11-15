//! Exact synthesis of fractional Gaussian noise (feature `synth`).
//!
//! Uses the Davies-Harte circulant embedding (Davies and Harte 1987): the
//! fGn autocovariance `g(k) = (|k+1|^{2H} - 2|k|^{2H} + |k-1|^{2H}) / 2` is
//! embedded in a circulant matrix of size `2n`, whose eigenvalues come from one
//! FFT. A second FFT of scaled complex Gaussian noise yields an exact sample.
//! Useful for testing estimators against a known Hurst exponent.

use rustfft::{num_complex::Complex, FftPlanner};

use crate::error::{Error, Result};

/// Small deterministic generator (SplitMix64 + Box-Muller) so samples are
/// reproducible from a seed without pulling in `rand`.
struct Gauss {
    state: u64,
    spare: Option<f64>,
}

impl Gauss {
    fn new(seed: u64) -> Self {
        Self {
            state: seed,
            spare: None,
        }
    }

    fn uniform(&mut self) -> f64 {
        self.state = self.state.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^= z >> 31;
        ((z >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }

    fn next(&mut self) -> f64 {
        if let Some(v) = self.spare.take() {
            return v;
        }
        let (u1, u2) = (self.uniform(), self.uniform());
        let r = (-2.0 * u1.ln()).sqrt();
        let t = 2.0 * core::f64::consts::PI * u2;
        self.spare = Some(r * t.sin());
        r * t.cos()
    }
}

/// `n` samples of unit-variance fractional Gaussian noise with Hurst
/// exponent `hurst` in `(0, 1)`, reproducible from `seed`.
pub fn fgn(n: usize, hurst: f64, seed: u64) -> Result<Vec<f64>> {
    if !(hurst > 0.0 && hurst < 1.0) {
        return Err(Error::InvalidParameter("hurst must be in (0, 1)"));
    }
    if n == 0 {
        return Ok(Vec::new());
    }
    let m = 2 * n;
    let h2 = 2.0 * hurst;
    let gamma = |k: f64| {
        0.5 * ((k + 1.0).abs().powf(h2) - 2.0 * k.abs().powf(h2) + (k - 1.0).abs().powf(h2))
    };
    let mut row: Vec<Complex<f64>> = (0..m)
        .map(|i| {
            let k = if i <= n { i } else { m - i };
            Complex::new(gamma(k as f64), 0.0)
        })
        .collect();
    let mut planner = FftPlanner::<f64>::new();
    let fft = planner.plan_fft_forward(m);
    fft.process(&mut row);
    let mut g = Gauss::new(seed);
    let mut w: Vec<Complex<f64>> = row
        .iter()
        .map(|lam| {
            let l = lam.re.max(0.0);
            let s = (l / m as f64).sqrt();
            Complex::new(s * g.next(), s * g.next())
        })
        .collect();
    fft.process(&mut w);
    Ok(w[..n].iter().map(|c| c.re).collect())
}

/// Fractional Brownian motion path: the cumulative sum of [`fgn`], starting
/// at the first increment.
pub fn fbm(n: usize, hurst: f64, seed: u64) -> Result<Vec<f64>> {
    let mut acc = 0.0;
    Ok(fgn(n, hurst, seed)?
        .into_iter()
        .map(|v| {
            acc += v;
            acc
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn variance_and_lag_one_correlation() {
        for &h in &[0.3, 0.5, 0.8] {
            let x = fgn(1 << 15, h, 42).unwrap();
            let n = x.len() as f64;
            let mean = x.iter().sum::<f64>() / n;
            let var = x.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n;
            let c1 = x
                .windows(2)
                .map(|w| (w[0] - mean) * (w[1] - mean))
                .sum::<f64>()
                / n
                / var;
            let want = 0.5 * (2f64.powf(2.0 * h) - 2.0);
            assert!((var - 1.0).abs() < 0.05, "H={h} var={var}");
            assert!((c1 - want).abs() < 0.03, "H={h} rho1={c1} want={want}");
        }
    }

    #[test]
    fn reproducible() {
        assert_eq!(fgn(64, 0.7, 9).unwrap(), fgn(64, 0.7, 9).unwrap());
        assert_ne!(fgn(64, 0.7, 9).unwrap(), fgn(64, 0.7, 10).unwrap());
        assert!(fgn(10, 1.0, 1).is_err());
    }
}
