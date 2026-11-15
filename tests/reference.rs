//! Numerical reference checks against independent Python implementations:
//! the `MFDFA` package (F_q(s)), `antropy.higuchi_fd` and PyWavelets
//! (`db4`, periodization). Expected values are produced by
//! `tools/gen_reference.py`; the input signals are rebuilt here bit for bit
//! from the same 64-bit LCG.

mod reference_data;
use reference_data::*;

use mfdfa::wavelet::dwt_level;
use mfdfa::{fluctuation_function, higuchi_fd, mfdfa, wavelet_hurst_detail};

const LCG_A: u64 = 6364136223846793005;
const LCG_C: u64 = 1442695040888963407;

fn lcg_uniform(n: usize, seed: u64) -> Vec<f64> {
    let mut s = seed;
    (0..n)
        .map(|_| {
            s = s.wrapping_mul(LCG_A).wrapping_add(LCG_C);
            (s >> 11) as f64 / (1u64 << 53) as f64
        })
        .collect()
}

fn cascade(levels: usize, p: f64, seed: u64) -> Vec<f64> {
    let mut s = seed;
    let mut v = vec![1.0];
    for _ in 0..levels {
        let mut next = Vec::with_capacity(2 * v.len());
        for x in v {
            s = s.wrapping_mul(LCG_A).wrapping_add(LCG_C);
            if (s >> 63) & 1 == 1 {
                next.push(x * p);
                next.push(x * (1.0 - p));
            } else {
                next.push(x * (1.0 - p));
                next.push(x * p);
            }
        }
        v = next;
    }
    v
}

fn signal_a() -> Vec<f64> {
    let n = 4096;
    let u = lcg_uniform(n, 12345);
    let c = cascade(12, 0.35, 777);
    u.iter()
        .zip(&c)
        .map(|(u, c)| (u - 0.5) * c * n as f64)
        .collect()
}

fn signal_b() -> Vec<f64> {
    lcg_uniform(3000, 2024).iter().map(|u| u - 0.5).collect()
}

fn cumsum(x: &[f64]) -> Vec<f64> {
    let mut acc = 0.0;
    x.iter()
        .map(|v| {
            acc += v;
            acc
        })
        .collect()
}

fn assert_close(got: &[f64], want: &[f64], rel: f64, what: &str) {
    assert_eq!(got.len(), want.len(), "{what}: length");
    for (i, (g, w)) in got.iter().zip(want).enumerate() {
        let tol = rel * w.abs().max(1.0);
        assert!(
            (g - w).abs() <= tol,
            "{what}[{i}]: got {g}, want {w} (diff {})",
            (g - w).abs()
        );
    }
}

#[test]
fn signals_are_rebuilt_exactly() {
    assert_eq!(&signal_a()[..8], &SIGNAL_A_HEAD);
    assert_eq!(&signal_b()[..8], &SIGNAL_B_HEAD);
}

#[allow(clippy::too_many_arguments)]
fn check_mfdfa(
    x: &[f64],
    order: usize,
    fq: &[f64],
    hq: &[f64],
    tau: &[f64],
    alpha: &[f64],
    fa: &[f64],
    tag: &str,
) {
    let got = fluctuation_function(x, &Q, &SCALES, order).unwrap();
    let flat: Vec<f64> = got.iter().flatten().copied().collect();
    assert_close(&flat, fq, 1e-10, &format!("{tag} F_q(s)"));
    let spec = mfdfa(x, &Q, &SCALES, order).unwrap();
    assert_close(&spec.h_q, hq, 1e-9, &format!("{tag} h(q)"));
    assert_close(&spec.tau_q, tau, 1e-9, &format!("{tag} tau(q)"));
    assert_close(&spec.alpha, alpha, 1e-9, &format!("{tag} alpha"));
    assert_close(&spec.f_alpha, fa, 1e-9, &format!("{tag} f(alpha)"));
}

#[test]
fn mfdfa_matches_python_mfdfa_package() {
    let (a, b) = (signal_a(), signal_b());
    check_mfdfa(
        &a,
        1,
        &FQ_A_O1,
        &HQ_A_O1,
        &TAU_A_O1,
        &ALPHA_A_O1,
        &FALPHA_A_O1,
        "A order 1",
    );
    check_mfdfa(
        &a,
        2,
        &FQ_A_O2,
        &HQ_A_O2,
        &TAU_A_O2,
        &ALPHA_A_O2,
        &FALPHA_A_O2,
        "A order 2",
    );
    check_mfdfa(
        &a,
        3,
        &FQ_A_O3,
        &HQ_A_O3,
        &TAU_A_O3,
        &ALPHA_A_O3,
        &FALPHA_A_O3,
        "A order 3",
    );
    check_mfdfa(
        &b,
        1,
        &FQ_B_O1,
        &HQ_B_O1,
        &TAU_B_O1,
        &ALPHA_B_O1,
        &FALPHA_B_O1,
        "B order 1",
    );
    check_mfdfa(
        &b,
        2,
        &FQ_B_O2,
        &HQ_B_O2,
        &TAU_B_O2,
        &ALPHA_B_O2,
        &FALPHA_B_O2,
        "B order 2",
    );
    check_mfdfa(
        &b,
        3,
        &FQ_B_O3,
        &HQ_B_O3,
        &TAU_B_O3,
        &ALPHA_B_O3,
        &FALPHA_B_O3,
        "B order 3",
    );
}

#[test]
fn cascade_signal_is_multifractal_and_noise_is_not() {
    let a = mfdfa(&signal_a(), &Q, &SCALES, 1).unwrap();
    let b = mfdfa(&signal_b(), &Q, &SCALES, 1).unwrap();
    assert!(
        a.width() > 3.0 * b.width(),
        "cascade width {} vs noise width {}",
        a.width(),
        b.width()
    );
    assert!((b.hurst() - 0.5).abs() < 0.06);
}

#[test]
fn higuchi_matches_antropy() {
    let (a, b) = (signal_a(), signal_b());
    let (pa, pb) = (cumsum(&a), cumsum(&b));
    for (i, &k) in KMAX.iter().enumerate() {
        assert_close(
            &[higuchi_fd(&pa, k)],
            &[HIGUCHI_PATH_A[i]],
            1e-10,
            &format!("path A k_max={k}"),
        );
        assert_close(
            &[higuchi_fd(&pb, k)],
            &[HIGUCHI_PATH_B[i]],
            1e-10,
            &format!("path B k_max={k}"),
        );
        assert_close(
            &[higuchi_fd(&b, k)],
            &[HIGUCHI_NOISE_B[i]],
            1e-10,
            &format!("noise B k_max={k}"),
        );
    }
}

#[test]
fn dwt_matches_pywavelets_periodization() {
    let a = signal_a();
    let (ca, cd) = dwt_level(&a[..256]);
    assert_close(&ca, &DWT_A256_APPROX, 1e-12, "db4 approximation");
    assert_close(&cd, &DWT_A256_DETAIL, 1e-12, "db4 detail");
}

#[test]
fn wavelet_hurst_matches_pywavelets_port() {
    for (x, h, ys, tag) in [
        (signal_a(), WAVELET_H_A[0], &WAVELET_LOG2E_A[..], "A"),
        (signal_b(), WAVELET_H_B[0], &WAVELET_LOG2E_B[..], "B"),
    ] {
        let w = wavelet_hurst_detail(&x, 1, None).unwrap();
        assert_close(&w.log2_energy, ys, 1e-10, &format!("{tag} log2 energy"));
        assert_close(&[w.hurst], &[h], 1e-10, &format!("{tag} wavelet H"));
    }
}
