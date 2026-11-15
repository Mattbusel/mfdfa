//! Orthonormal Daubechies-4 (`db4`, 8 taps) discrete wavelet transform and
//! the wavelet (Abry-Veitch) Hurst estimator.
//!
//! The transform uses periodic boundaries with the same filter alignment as
//! PyWavelets' `mode="periodization"`, so for a power-of-two input
//! [`dwt_level`] returns exactly `pywt.dwt(x, "db4", mode="periodization")`
//! (checked in `tests/reference.rs`). Because the periodic transform is
//! orthogonal, [`inverse_dwt`] is its transpose and reconstructs the input to
//! rounding error.
//!
//! # Wavelet Hurst
//!
//! For a signal with stationary increments of Hurst exponent `H` (fractional
//! Gaussian noise, or log returns with long memory) the detail coefficient
//! variance at octave `j` scales as `2^{j (2H - 1)}`. [`wavelet_hurst`] fits
//! `log2 mean(d_j^2)` against `j` by weighted least squares (weights = number
//! of coefficients per octave, Abry and Veitch 1998) and returns
//! `H = (slope + 1) / 2`. Pass increments (returns), not the price path.

use crate::error::{check_finite, wls_slope, Error, Result};

/// db4 decomposition low-pass filter (PyWavelets `dec_lo`).
pub const DB4_DEC_LO: [f64; 8] = [
    -0.010597401785069032,
    0.0328830116668852,
    0.030841381835560764,
    -0.18703481171909309,
    -0.027983769416859854,
    0.6308807679298589,
    0.7148465705529157,
    0.2303778133088965,
];

/// db4 decomposition high-pass filter (PyWavelets `dec_hi`).
pub const DB4_DEC_HI: [f64; 8] = [
    -0.2303778133088965,
    0.7148465705529157,
    -0.6308807679298589,
    -0.027983769416859854,
    0.18703481171909309,
    0.030841381835560764,
    -0.0328830116668852,
    -0.010597401785069032,
];

const SHIFT: usize = DB4_DEC_LO.len() / 2;

/// Result of a multi-level DWT.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct WaveletDecomp {
    /// Final coarse approximation (scaling) coefficients.
    pub approximation: Vec<f64>,
    /// Detail (wavelet) coefficients per level, finest first.
    pub details: Vec<Vec<f64>>,
    /// Original signal length (before zero padding).
    pub original_len: usize,
    /// Number of decomposition levels.
    pub levels: usize,
}

impl WaveletDecomp {
    /// Per-scale energy: `[approximation, detail_0 (finest), detail_1, ...]`.
    pub fn energy(&self) -> Vec<f64> {
        let mut e = Vec::with_capacity(self.levels + 1);
        e.push(l2_energy(&self.approximation));
        for d in &self.details {
            e.push(l2_energy(d));
        }
        e
    }

    /// Fraction of total energy at each scale (same order as [`Self::energy`]).
    pub fn energy_fractions(&self) -> Vec<f64> {
        let e = self.energy();
        let total: f64 = e.iter().sum();
        if total < 1e-16 {
            return vec![0.0; e.len()];
        }
        e.iter().map(|v| v / total).collect()
    }

    /// Index of the detail level with the most energy (0 = finest).
    pub fn dominant_detail_scale(&self) -> usize {
        self.details
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| {
                l2_energy(a)
                    .partial_cmp(&l2_energy(b))
                    .unwrap_or(core::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0)
    }
}

fn l2_energy(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum()
}

/// One level of the periodic db4 DWT. `signal.len()` must be even and at
/// least 2. Returns `(approximation, detail)`, each of half the length.
///
/// Coefficient `i` is `sum_k h[k] x[(2 i + 4 - k) mod n]`, the PyWavelets
/// periodization convention.
pub fn dwt_level(signal: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let n = signal.len();
    assert!(n >= 2 && n % 2 == 0, "dwt_level needs an even length >= 2");
    let half = n / 2;
    let mut a = vec![0.0; half];
    let mut d = vec![0.0; half];
    for i in 0..half {
        let mut lo = 0.0;
        let mut hi = 0.0;
        for k in 0..DB4_DEC_LO.len() {
            let idx = (2 * i + SHIFT + n * DB4_DEC_LO.len() - k) % n;
            lo += DB4_DEC_LO[k] * signal[idx];
            hi += DB4_DEC_HI[k] * signal[idx];
        }
        a[i] = lo;
        d[i] = hi;
    }
    (a, d)
}

/// Exact inverse of [`dwt_level`].
pub fn idwt_level(approx: &[f64], detail: &[f64]) -> Vec<f64> {
    assert_eq!(approx.len(), detail.len());
    let half = approx.len();
    let n = 2 * half;
    let mut out = vec![0.0; n];
    for i in 0..half {
        for k in 0..DB4_DEC_LO.len() {
            let idx = (2 * i + SHIFT + n * DB4_DEC_LO.len() - k) % n;
            out[idx] += DB4_DEC_LO[k] * approx[i] + DB4_DEC_HI[k] * detail[i];
        }
    }
    out
}

/// Multi-level DWT. The signal is zero padded to the next power of two (at
/// least 8) and `levels` is capped so the coarsest level keeps at least 8
/// coefficients.
pub fn forward_dwt(signal: &[f64], levels: usize) -> WaveletDecomp {
    let original_len = signal.len();
    let padded_len = original_len.max(8).next_power_of_two();
    let mut current = signal.to_vec();
    current.resize(padded_len, 0.0);
    let levels = levels.min(max_levels(padded_len));
    let mut details = Vec::with_capacity(levels);
    for _ in 0..levels {
        let (a, d) = dwt_level(&current);
        details.push(d);
        current = a;
    }
    WaveletDecomp {
        approximation: current,
        details,
        original_len,
        levels,
    }
}

/// Reconstruct the signal from [`forward_dwt`] output (trimmed to the
/// original length).
pub fn inverse_dwt(decomp: &WaveletDecomp) -> Vec<f64> {
    let mut current = decomp.approximation.clone();
    for d in decomp.details.iter().rev() {
        current = idwt_level(&current, d);
    }
    current.truncate(decomp.original_len);
    current
}

fn max_levels(n: usize) -> usize {
    let mut levels = 0;
    let mut size = n;
    while size >= 2 * DB4_DEC_LO.len() {
        size /= 2;
        levels += 1;
    }
    levels.max(1)
}

/// Detailed output of [`wavelet_hurst_detail`].
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct WaveletHurst {
    /// Estimated Hurst exponent `(slope + 1) / 2`.
    pub hurst: f64,
    /// Octaves used in the fit (1 = finest).
    pub octaves: Vec<usize>,
    /// `log2` of the mean squared detail coefficient per octave.
    pub log2_energy: Vec<f64>,
    /// Number of detail coefficients per octave (the regression weights).
    pub counts: Vec<usize>,
}

/// Wavelet Hurst estimate over octaves `min_octave..=max_octave`.
///
/// `max_octave = None` uses every octave that still has at least 8
/// coefficients. At each level an odd-length approximation drops its last
/// sample, so any length works. Needs at least two octaves.
pub fn wavelet_hurst_detail(
    x: &[f64],
    min_octave: usize,
    max_octave: Option<usize>,
) -> Result<WaveletHurst> {
    check_finite(x)?;
    if min_octave == 0 {
        return Err(Error::InvalidParameter("octaves start at 1"));
    }
    let mut current: Vec<f64> = x.to_vec();
    // Detail energy below this is rounding noise (constant or polynomial input).
    let floor = 1e-24 * x.iter().map(|v| v * v).sum::<f64>() / x.len().max(1) as f64;
    let mut octaves = Vec::new();
    let mut log2_energy = Vec::new();
    let mut counts = Vec::new();
    let mut j = 0;
    loop {
        if current.len() % 2 == 1 {
            current.pop();
        }
        if current.len() < 16 {
            break;
        }
        j += 1;
        if let Some(mx) = max_octave {
            if j > mx {
                break;
            }
        }
        let (a, d) = dwt_level(&current);
        if j >= min_octave {
            let e = l2_energy(&d) / d.len() as f64;
            if e.is_nan() || e <= floor || e == 0.0 {
                return Err(Error::Degenerate("zero detail energy (constant series)"));
            }
            octaves.push(j);
            log2_energy.push(e.log2());
            counts.push(d.len());
        }
        current = a;
    }
    if octaves.len() < 2 {
        let needed = 16usize << min_octave.max(1);
        return Err(Error::TooShort {
            needed,
            got: x.len(),
        });
    }
    let xs: Vec<f64> = octaves.iter().map(|&j| j as f64).collect();
    let w: Vec<f64> = counts.iter().map(|&c| c as f64).collect();
    let slope = wls_slope(&xs, &log2_energy, &w);
    Ok(WaveletHurst {
        hurst: (slope + 1.0) / 2.0,
        octaves,
        log2_energy,
        counts,
    })
}

/// Wavelet (Abry-Veitch) Hurst exponent of an increment series, using every
/// octave from the finest down to the last one with at least 8 coefficients.
///
/// ```
/// let mut s = 1u64;
/// let noise: Vec<f64> = (0..4096)
///     .map(|_| {
///         s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
///         (s >> 11) as f64 / (1u64 << 53) as f64 - 0.5
///     })
///     .collect();
/// let h = mfdfa::wavelet_hurst(&noise).unwrap();
/// assert!((h - 0.5).abs() < 0.08);
/// ```
pub fn wavelet_hurst(x: &[f64]) -> Result<f64> {
    wavelet_hurst_detail(x, 1, None).map(|w| w.hurst)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sine_signal(n: usize, freq: f64) -> Vec<f64> {
        (0..n)
            .map(|i| (2.0 * core::f64::consts::PI * freq * i as f64 / n as f64).sin())
            .collect()
    }

    #[test]
    fn filters_are_orthonormal() {
        let lo: f64 = DB4_DEC_LO.iter().map(|v| v * v).sum();
        let hi: f64 = DB4_DEC_HI.iter().map(|v| v * v).sum();
        let cross: f64 = DB4_DEC_LO.iter().zip(&DB4_DEC_HI).map(|(a, b)| a * b).sum();
        assert!((lo - 1.0).abs() < 1e-12 && (hi - 1.0).abs() < 1e-12 && cross.abs() < 1e-12);
        assert!((DB4_DEC_LO.iter().sum::<f64>() - 2f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn forward_dwt_returns_correct_structure() {
        let decomp = forward_dwt(&sine_signal(128, 5.0), 3);
        assert_eq!(decomp.levels, 3);
        assert_eq!(decomp.details.len(), 3);
        assert_eq!(decomp.approximation.len(), 16);
    }

    #[test]
    fn energy_fractions_sum_to_one() {
        let fracs = forward_dwt(&sine_signal(128, 3.0), 4).energy_fractions();
        assert!((fracs.iter().sum::<f64>() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn parseval_holds() {
        let signal = sine_signal(256, 7.0);
        let e_in: f64 = l2_energy(&signal);
        let e_out: f64 = forward_dwt(&signal, 4).energy().iter().sum();
        assert!((e_in - e_out).abs() < 1e-9 * e_in);
    }

    #[test]
    fn perfect_reconstruction() {
        for n in [64usize, 100, 257] {
            let signal: Vec<f64> = (0..n)
                .map(|i| ((i * 37 % 11) as f64).sin() + i as f64 * 0.01)
                .collect();
            let rec = inverse_dwt(&forward_dwt(&signal, 5));
            assert_eq!(rec.len(), n);
            for (a, b) in signal.iter().zip(&rec) {
                assert!((a - b).abs() < 1e-10, "n={n}: {a} vs {b}");
            }
        }
    }

    #[test]
    fn dominant_scale_is_valid() {
        let decomp = forward_dwt(&sine_signal(128, 8.0), 4);
        assert!(decomp.dominant_detail_scale() < decomp.levels);
    }

    #[test]
    fn wavelet_hurst_errors() {
        let short: Vec<f64> = (0..20).map(|i| ((i * i) % 7) as f64).collect();
        assert!(matches!(wavelet_hurst(&short), Err(Error::TooShort { .. })));
        assert!(matches!(
            wavelet_hurst(&[1.0; 256]),
            Err(Error::Degenerate(_))
        ));
        assert!(matches!(
            wavelet_hurst_detail(&[1.0; 256], 0, None),
            Err(Error::InvalidParameter(_))
        ));
    }
}
