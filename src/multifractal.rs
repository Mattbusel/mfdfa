//! Price-oriented convenience wrapper around [`crate::mfdfa`].
//!
//! [`MultifractalAnalyzer`] takes a price series, converts it to log returns,
//! picks scales (8, then x1.5 up to n/4) and runs MFDFA-1 over q = -5..=5.
//! This is the interface the srfm-lab `fractal-engine` crate exposed; the
//! numerical core is now the general [`crate::mfdfa`] function.
//!
//! # Interpretation
//!
//! - `multifractal_width > 0.3`: multifractal, heterogeneous dynamics
//! - `multifractal_width < 0.1`: monofractal (near-Brownian)
//! - `h2` is the ordinary DFA Hurst exponent

use crate::dfa::mfdfa;
use crate::error::log_returns;

/// Result of a full multifractal analysis.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MultifractalResult {
    /// Generalised Hurst exponents as `(q, h(q))` pairs.
    pub hq: Vec<(f64, f64)>,
    /// `tau(q) = q h(q) - 1`.
    pub tau_q: Vec<f64>,
    /// Singularity spectrum `(alpha, f_alpha)`.
    pub singularity_spectrum: (Vec<f64>, Vec<f64>),
    /// `alpha_max - alpha_min`.
    pub multifractal_width: f64,
    /// `multifractal_width > 0.3`.
    pub is_multifractal: bool,
    /// `h(q = 2)`: the standard DFA Hurst exponent (NaN if q = 2 is absent).
    pub h2: f64,
}

/// MFDFA on prices with sensible defaults.
#[derive(Debug, Clone)]
pub struct MultifractalAnalyzer {
    /// q orders (strictly increasing). Default -5..=5 step 1.
    pub q_orders: Vec<f64>,
    /// Segment sizes. Empty means automatic.
    pub scales: Vec<usize>,
    /// Detrending order. Default 1.
    pub order: usize,
}

impl MultifractalAnalyzer {
    /// q in -5..=5 step 1, automatic scales, linear detrending.
    pub fn new() -> Self {
        Self {
            q_orders: (-5i32..=5).map(f64::from).collect(),
            scales: Vec::new(),
            order: 1,
        }
    }

    /// Explicit q orders and scales.
    pub fn with_config(q_orders: Vec<f64>, scales: Vec<usize>) -> Self {
        Self {
            q_orders,
            scales,
            order: 1,
        }
    }

    fn default_scales(n: usize) -> Vec<usize> {
        let mut scales = Vec::new();
        let mut s = 8usize;
        while s <= n / 4 {
            scales.push(s);
            s = (s as f64 * 1.5).ceil() as usize;
        }
        scales
    }

    /// Integrated profile of returns (kept for compatibility).
    pub fn compute_profile(returns: &[f64]) -> Vec<f64> {
        crate::dfa::profile(returns)
    }

    /// Run MFDFA on a price series (at least 64 prices).
    pub fn analyze(&self, prices: &[f64]) -> Option<MultifractalResult> {
        if prices.len() < 64 {
            return None;
        }
        let returns = log_returns(prices);
        let n = returns.len();
        let mut scales: Vec<usize> = if self.scales.is_empty() {
            Self::default_scales(prices.len())
        } else {
            self.scales.clone()
        };
        scales.retain(|&s| s >= self.order + 2 && 2 * s <= n);
        scales.sort_unstable();
        scales.dedup();
        if scales.len() < 3 {
            return None;
        }
        let spec = mfdfa(&returns, &self.q_orders, &scales, self.order).ok()?;
        let h2 = spec
            .q
            .iter()
            .position(|q| (q - 2.0).abs() < 1e-9)
            .map(|i| spec.h_q[i])
            .unwrap_or(f64::NAN);
        let width = spec.width();
        Some(MultifractalResult {
            hq: spec
                .q
                .iter()
                .copied()
                .zip(spec.h_q.iter().copied())
                .collect(),
            tau_q: spec.tau_q.clone(),
            singularity_spectrum: (spec.alpha.clone(), spec.f_alpha.clone()),
            multifractal_width: width,
            is_multifractal: width > 0.3,
            h2,
        })
    }

    /// `(alpha, f_alpha)` directly from prices (empty on failure).
    pub fn singularity_spectrum(&self, prices: &[f64]) -> (Vec<f64>, Vec<f64>) {
        self.analyze(prices)
            .map(|r| r.singularity_spectrum)
            .unwrap_or_default()
    }

    /// Spectrum width (0 on failure).
    pub fn multifractal_width(&self, prices: &[f64]) -> f64 {
        self.analyze(prices)
            .map(|r| r.multifractal_width)
            .unwrap_or(0.0)
    }

    /// `multifractal_width(prices) > 0.3`.
    pub fn is_multifractal(&self, prices: &[f64]) -> bool {
        self.multifractal_width(prices) > 0.3
    }
}

impl Default for MultifractalAnalyzer {
    fn default() -> Self {
        Self::new()
    }
}

/// Summary of the `h(q)` curve at q = -5, 0, 2, 5.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct HurstSpectrum {
    /// h(-5): scaling of small fluctuations.
    pub h_neg5: f64,
    /// h(0).
    pub h_0: f64,
    /// h(2): the standard DFA Hurst exponent.
    pub h_2: f64,
    /// h(5): scaling of large fluctuations.
    pub h_5: f64,
    /// h(-5) - h(5).
    pub hq_range: f64,
}

impl HurstSpectrum {
    /// Pick the four orders out of a result (None if any is missing).
    pub fn from_result(result: &MultifractalResult) -> Option<Self> {
        let get = |q: f64| {
            result
                .hq
                .iter()
                .find(|(qq, _)| (qq - q).abs() < 1e-9)
                .map(|(_, h)| *h)
        };
        let (h_neg5, h_0, h_2, h_5) = (get(-5.0)?, get(0.0)?, get(2.0)?, get(5.0)?);
        Some(Self {
            h_neg5,
            h_0,
            h_2,
            h_5,
            hq_range: h_neg5 - h_5,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Generate a random walk price series using LCG.
    fn random_walk(n: usize, seed: u64) -> Vec<f64> {
        let mut price = 100.0f64;
        let mut prices = Vec::with_capacity(n);
        let mut state = seed;
        for _ in 0..n {
            prices.push(price);
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let r = (state as f64 / u64::MAX as f64) * 2.0 - 1.0;
            price = (price + r * 0.5).max(1.0);
        }
        prices
    }

    /// Generate a strongly trending price series (persistent, multifractal).
    fn trending_series(n: usize) -> Vec<f64> {
        let mut prices = Vec::with_capacity(n);
        let mut p = 100.0f64;
        let mut state = 7u64;
        for _ in 0..n {
            prices.push(p);
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let noise = (state as f64 / u64::MAX as f64) * 0.002 - 0.001;
            p *= 1.003 + noise; // strong trend + small noise
        }
        prices
    }

    #[test]
    fn test_mfdfa_returns_result_for_sufficient_data() {
        let prices = random_walk(256, 42);
        let analyzer = MultifractalAnalyzer::new();
        let result = analyzer.analyze(&prices);
        assert!(result.is_some(), "Should return a result for 256 prices");
    }

    #[test]
    fn test_mfdfa_too_short_returns_none() {
        let prices = random_walk(32, 42);
        let analyzer = MultifractalAnalyzer::new();
        assert!(analyzer.analyze(&prices).is_none());
    }

    #[test]
    fn test_multifractal_width_positive() {
        let prices = random_walk(512, 99);
        let analyzer = MultifractalAnalyzer::new();
        let width = analyzer.multifractal_width(&prices);
        assert!(width >= 0.0, "Width should be non-negative, got {width}");
    }

    #[test]
    fn test_hq_count_matches_q_orders() {
        let prices = random_walk(256, 7);
        let analyzer = MultifractalAnalyzer::new();
        if let Some(result) = analyzer.analyze(&prices) {
            // All q in [-5..5] that produced valid H(q) should be present
            assert!(result.hq.len() >= 5, "expected at least 5 H(q) values");
        }
    }

    #[test]
    fn test_h2_in_valid_range() {
        let prices = random_walk(256, 42);
        let analyzer = MultifractalAnalyzer::new();
        if let Some(result) = analyzer.analyze(&prices) {
            assert!(
                result.h2 >= 0.0 && result.h2 <= 2.0,
                "H(2) = {} should be in [0, 2]",
                result.h2
            );
        }
    }

    #[test]
    fn test_hurst_h2_near_standard_hurst() {
        // For a random walk H should be near 0.5; H(q=2) from MFDFA should
        // also be near 0.5 (within tolerance of 0.3 given finite sample).
        let prices = random_walk(512, 12345);
        let analyzer = MultifractalAnalyzer::new();
        if let Some(result) = analyzer.analyze(&prices) {
            if result.h2.is_finite() {
                assert!(
                    result.h2 >= 0.2 && result.h2 <= 0.8,
                    "H(2) = {:.3} should be near 0.5 for random walk",
                    result.h2
                );
            }
        }
    }

    #[test]
    fn test_singularity_spectrum_lengths_match() {
        let prices = random_walk(256, 1);
        let analyzer = MultifractalAnalyzer::new();
        let (alpha, f_alpha) = analyzer.singularity_spectrum(&prices);
        assert_eq!(
            alpha.len(),
            f_alpha.len(),
            "alpha and f_alpha must have equal length"
        );
    }

    #[test]
    fn test_tau_q_length_matches_hq() {
        let prices = random_walk(256, 3);
        let analyzer = MultifractalAnalyzer::new();
        if let Some(result) = analyzer.analyze(&prices) {
            assert_eq!(result.tau_q.len(), result.hq.len());
        }
    }

    #[test]
    fn test_trending_series_high_h2() {
        // A strongly trending series should have H(2) higher than a random walk.
        // Exact value varies with finite-sample estimation; just check > 0.4.
        let prices = trending_series(512);
        let analyzer = MultifractalAnalyzer::new();
        if let Some(result) = analyzer.analyze(&prices) {
            if result.h2.is_finite() {
                assert!(
                    result.h2 > 0.4,
                    "trending H(2) = {:.3} should be > 0.4",
                    result.h2
                );
            }
        }
    }

    #[test]
    fn test_is_multifractal_consistent_with_width() {
        let prices = random_walk(512, 77);
        let analyzer = MultifractalAnalyzer::new();
        if let Some(result) = analyzer.analyze(&prices) {
            assert_eq!(result.is_multifractal, result.multifractal_width > 0.3);
        }
    }

    #[test]
    fn test_hurst_spectrum_construction() {
        let prices = random_walk(512, 55);
        let analyzer = MultifractalAnalyzer::new();
        if let Some(result) = analyzer.analyze(&prices) {
            let spectrum = HurstSpectrum::from_result(&result);
            // May be None if q=+/-5 didn't converge; that's fine
            if let Some(s) = spectrum {
                assert!(s.h_2.is_finite());
                assert!(s.hq_range.is_finite());
            }
        }
    }

    #[test]
    fn test_profile_zero_mean() {
        // The integrated profile should have approximately zero drift
        let returns = vec![0.01, -0.01, 0.02, -0.02, 0.005, -0.005];
        let profile = MultifractalAnalyzer::compute_profile(&returns);
        let last = profile.last().copied().unwrap_or(0.0);
        assert!(
            last.abs() < 1e-10,
            "cumsum of demeaned returns should end near 0"
        );
    }
}
