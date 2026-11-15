//! Rolling, bar-by-bar estimators.
//!
//! [`RollingHurst`] and [`RollingHiguchi`] keep a fixed window and return a
//! fresh estimate on every update (NaN until the window is full or when the
//! window is degenerate, for example flat prices). With the `ta` feature they
//! implement `ta::Next<f64>`, `ta::Next<&T: ta::Close>`, `ta::Reset` and
//! `ta::Period`, so they drop into a `ta` indicator pipeline.

use std::collections::VecDeque;
use std::fmt;

use crate::dfa::{default_scales, dfa_hurst};
use crate::error::{Error, Result};
use crate::higuchi::higuchi_fd;

/// What a [`RollingHurst`] is fed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Input {
    /// Positive prices; log returns are taken internally.
    Prices,
    /// Increments directly (returns, order flow, any stationary signal).
    Increments,
}

/// Rolling DFA Hurst exponent over the last `period` increments.
#[derive(Debug, Clone)]
pub struct RollingHurst {
    period: usize,
    order: usize,
    input: Input,
    scales: Vec<usize>,
    last_price: Option<f64>,
    buf: VecDeque<f64>,
}

impl RollingHurst {
    /// Rolling Hurst on prices (log returns of the last `period + 1` prices),
    /// DFA-1. `period` must be at least 48.
    pub fn new(period: usize) -> Result<Self> {
        Self::with_options(period, 1, Input::Prices)
    }

    /// Rolling Hurst fed increments directly, DFA-1.
    pub fn on_increments(period: usize) -> Result<Self> {
        Self::with_options(period, 1, Input::Increments)
    }

    /// Full control: detrending `order` and input kind.
    pub fn with_options(period: usize, order: usize, input: Input) -> Result<Self> {
        let scales = default_scales(period, order);
        if scales.len() < 3 {
            return Err(Error::InvalidParameter(
                "period too short for at least three DFA scales (use 64 or more)",
            ));
        }
        Ok(Self {
            period,
            order,
            input,
            scales,
            last_price: None,
            buf: VecDeque::with_capacity(period + 1),
        })
    }

    /// Window length in increments.
    pub fn period(&self) -> usize {
        self.period
    }

    /// Push one observation and return the current estimate (NaN while
    /// warming up).
    pub fn update(&mut self, value: f64) -> f64 {
        let inc = match self.input {
            Input::Increments => Some(value),
            Input::Prices => {
                let r = match self.last_price {
                    Some(p) if p > 0.0 && value > 0.0 => Some((value / p).ln()),
                    Some(_) => Some(0.0),
                    None => None,
                };
                self.last_price = Some(value);
                r
            }
        };
        if let Some(v) = inc {
            if self.buf.len() == self.period {
                self.buf.pop_front();
            }
            self.buf.push_back(v);
        }
        self.value()
    }

    /// The estimate for the current window, without pushing anything.
    pub fn value(&self) -> f64 {
        if self.buf.len() < self.period {
            return f64::NAN;
        }
        let (a, b) = self.buf.as_slices();
        let w: Vec<f64> = a.iter().chain(b).copied().collect();
        dfa_hurst(&w, &self.scales, self.order).unwrap_or(f64::NAN)
    }

    /// Clear the window.
    pub fn clear(&mut self) {
        self.buf.clear();
        self.last_price = None;
    }
}

impl fmt::Display for RollingHurst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "HURST({})", self.period)
    }
}

/// Rolling Higuchi fractal dimension over the last `period` values of a curve
/// (prices or any path).
#[derive(Debug, Clone)]
pub struct RollingHiguchi {
    period: usize,
    k_max: usize,
    buf: VecDeque<f64>,
}

impl RollingHiguchi {
    /// `period` must be at least `2 * k_max` and `k_max` at least 2.
    pub fn new(period: usize, k_max: usize) -> Result<Self> {
        if k_max < 2 {
            return Err(Error::InvalidParameter("k_max must be at least 2"));
        }
        if period < 2 * k_max {
            return Err(Error::InvalidParameter("period must be at least 2 * k_max"));
        }
        Ok(Self {
            period,
            k_max,
            buf: VecDeque::with_capacity(period),
        })
    }

    /// Window length.
    pub fn period(&self) -> usize {
        self.period
    }

    /// Push one value and return the current estimate (NaN while warming up).
    pub fn update(&mut self, value: f64) -> f64 {
        if self.buf.len() == self.period {
            self.buf.pop_front();
        }
        self.buf.push_back(value);
        self.value()
    }

    /// The estimate for the current window.
    pub fn value(&self) -> f64 {
        if self.buf.len() < self.period {
            return f64::NAN;
        }
        let (a, b) = self.buf.as_slices();
        let w: Vec<f64> = a.iter().chain(b).copied().collect();
        higuchi_fd(&w, self.k_max)
    }

    /// Clear the window.
    pub fn clear(&mut self) {
        self.buf.clear();
    }
}

impl fmt::Display for RollingHiguchi {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "HIGUCHI({}, {})", self.period, self.k_max)
    }
}

#[cfg(feature = "ta")]
mod ta_impls {
    use super::{RollingHiguchi, RollingHurst};
    use ta::{Close, Next, Period, Reset};

    impl Next<f64> for RollingHurst {
        type Output = f64;
        fn next(&mut self, input: f64) -> f64 {
            self.update(input)
        }
    }

    impl<T: Close> Next<&T> for RollingHurst {
        type Output = f64;
        fn next(&mut self, input: &T) -> f64 {
            self.update(input.close())
        }
    }

    impl Reset for RollingHurst {
        fn reset(&mut self) {
            self.clear();
        }
    }

    impl Period for RollingHurst {
        fn period(&self) -> usize {
            self.period
        }
    }

    impl Next<f64> for RollingHiguchi {
        type Output = f64;
        fn next(&mut self, input: f64) -> f64 {
            self.update(input)
        }
    }

    impl<T: Close> Next<&T> for RollingHiguchi {
        type Output = f64;
        fn next(&mut self, input: &T) -> f64 {
            self.update(input.close())
        }
    }

    impl Reset for RollingHiguchi {
        fn reset(&mut self) {
            self.clear();
        }
    }

    impl Period for RollingHiguchi {
        fn period(&self) -> usize {
            self.period
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::higuchi::higuchi_fd;

    fn walk(n: usize) -> Vec<f64> {
        let mut s = 5u64;
        let mut p = 100.0;
        (0..n)
            .map(|_| {
                s = s
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                p *= 1.0 + 0.01 * ((s >> 11) as f64 / (1u64 << 53) as f64 - 0.5);
                p
            })
            .collect()
    }

    #[test]
    fn rolling_hurst_matches_batch() {
        let prices = walk(300);
        let mut ind = RollingHurst::new(128).unwrap();
        let out: Vec<f64> = prices.iter().map(|&p| ind.update(p)).collect();
        assert!(out[..128].iter().all(|v| v.is_nan()));
        let rets = crate::error::log_returns(&prices);
        let batch = dfa_hurst(&rets[rets.len() - 128..], &default_scales(128, 1), 1).unwrap();
        assert!((out[299] - batch).abs() < 1e-12);
        ind.clear();
        assert!(ind.update(100.0).is_nan());
    }

    #[test]
    fn rolling_higuchi_matches_batch() {
        let prices = walk(200);
        let mut ind = RollingHiguchi::new(64, 8).unwrap();
        let out: Vec<f64> = prices.iter().map(|&p| ind.update(p)).collect();
        assert!(out[..63].iter().all(|v| v.is_nan()));
        assert!((out[199] - higuchi_fd(&prices[136..], 8)).abs() < 1e-12);
    }

    #[test]
    fn constructor_validation() {
        assert!(RollingHurst::new(20).is_err());
        assert!(RollingHiguchi::new(10, 8).is_err());
        assert!(RollingHiguchi::new(64, 1).is_err());
        assert_eq!(RollingHurst::new(64).unwrap().to_string(), "HURST(64)");
    }
}
