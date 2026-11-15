//! Estimators against signals whose answer is known in closed form.

use mfdfa::{default_scales, dfa_hurst, higuchi_fd, mfdfa, wavelet_hurst};

/// Binomial p-model cascade (Kantelhardt et al. 2002, section 4).
fn binomial_cascade(levels: usize, p: f64) -> Vec<f64> {
    let mut v = vec![1.0];
    for _ in 0..levels {
        v = v.iter().flat_map(|x| [x * p, x * (1.0 - p)]).collect();
    }
    v
}

/// h(q) = 1/q - ln(p^q + (1-p)^q) / (q ln 2).
fn cascade_h(q: f64, p: f64) -> f64 {
    1.0 / q - (p.powf(q) + (1.0 - p).powf(q)).ln() / (q * 2f64.ln())
}

#[test]
fn binomial_cascade_matches_analytic_h_of_q() {
    // The shape of h(q) must match the closed form for every detrending
    // order; the level carries a small finite-size offset (about 0.03 here).
    let p = 0.3;
    let x = binomial_cascade(16, p);
    let q: Vec<f64> = vec![-3.0, -2.0, -1.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    let scales: Vec<usize> = (5..=12).map(|e| 1usize << e).collect();
    let i2 = q.iter().position(|v| *v == 2.0).unwrap();
    for order in 1..=3 {
        let spec = mfdfa(&x, &q, &scales, order).unwrap();
        for (i, qi) in q.iter().enumerate() {
            let want = cascade_h(*qi, p);
            let got = spec.h_q[i];
            let shape = (got - spec.h_q[i2]) - (want - cascade_h(2.0, p));
            assert!(
                shape.abs() < 0.005,
                "order {order}, q={qi}: shape error {shape:.5}"
            );
            assert!(
                (got - want).abs() < 0.08,
                "order {order}, q={qi}: h={got:.4}, analytic {want:.4}"
            );
        }
        // Spectrum width of the p-model: alpha_max - alpha_min = |log2((1-p)/p)| in the limit.
        assert!(spec.width() > 0.8, "order {order}: width {}", spec.width());
    }
}

#[test]
fn straight_line_has_dimension_one() {
    let x: Vec<f64> = (0..1000).map(|i| 0.5 * i as f64).collect();
    assert!((higuchi_fd(&x, 16) - 1.0).abs() < 1e-9);
}

#[cfg(feature = "synth")]
mod fgn {
    use super::*;
    use mfdfa::synth::{fbm, fgn};

    #[test]
    fn dfa_recovers_hurst() {
        for &h in &[0.3, 0.5, 0.7, 0.85] {
            let mut est = 0.0;
            for seed in 0..4 {
                let x = fgn(1 << 14, h, seed).unwrap();
                est += dfa_hurst(&x, &default_scales(x.len(), 1), 1).unwrap() / 4.0;
            }
            assert!((est - h).abs() < 0.05, "H={h}: DFA estimate {est:.3}");
        }
    }

    #[test]
    fn wavelet_recovers_hurst() {
        for &h in &[0.3, 0.5, 0.7, 0.85] {
            let mut est = 0.0;
            for seed in 0..4 {
                est += wavelet_hurst(&fgn(1 << 14, h, 100 + seed).unwrap()).unwrap() / 4.0;
            }
            assert!((est - h).abs() < 0.05, "H={h}: wavelet estimate {est:.3}");
        }
    }

    #[test]
    fn higuchi_on_fbm_is_two_minus_h() {
        for &h in &[0.3, 0.5, 0.7, 0.85] {
            let mut est = 0.0;
            for seed in 0..4 {
                est += higuchi_fd(&fbm(1 << 14, h, 200 + seed).unwrap(), 10) / 4.0;
            }
            assert!(
                (est - (2.0 - h)).abs() < 0.05,
                "H={h}: Higuchi {est:.3}, want {:.3}",
                2.0 - h
            );
        }
    }

    #[test]
    fn fgn_is_monofractal() {
        let x = fgn(1 << 14, 0.7, 7).unwrap();
        let spec = mfdfa(&x, &mfdfa::default_q(), &default_scales(x.len(), 1), 1).unwrap();
        assert!((spec.hurst() - 0.7).abs() < 0.06, "h(2) = {}", spec.hurst());
        assert!(spec.delta_h().abs() < 0.25, "delta h = {}", spec.delta_h());
    }
}
