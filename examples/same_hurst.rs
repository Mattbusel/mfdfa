//! Two return series with the same Hurst exponent and very different risk.
//!
//! cargo run --example same_hurst
//!
//! Both series are uncorrelated returns, so DFA gives both a Hurst exponent
//! near 0.5. The second one gets its volatility from a multiplicative
//! cascade (Mandelbrot's binomial measure, p = 0.25): calm stretches and
//! bursts nested inside bursts. A single Hurst number cannot tell them apart.
//! The width of the MFDFA spectrum can.

use mfdfa::{default_q, default_scales, mfdfa};

fn main() {
    let mut s = 7u64;
    let mut gauss = move || {
        // Box-Muller on a 64-bit LCG: reproducible, no dependencies.
        let mut u = || {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((s >> 11) as f64 + 0.5) / (1u64 << 53) as f64
        };
        let (a, b) = (u(), u());
        (-2.0 * a.ln()).sqrt() * (std::f64::consts::TAU * b).cos()
    };

    let levels = 13; // 2^13 = 8192 returns
    let n = 1usize << levels;
    let noise: Vec<f64> = (0..n).map(|_| gauss()).collect();

    // Binomial cascade: each interval hands share p to one half and 1 - p to
    // the other, which half chosen by the bits of the index.
    let p = 0.25;
    let weight: Vec<f64> = (0..n)
        .map(|i| {
            (0..levels)
                .map(|bit| {
                    if (i >> bit) & 1 == 0 {
                        2.0 * p
                    } else {
                        2.0 * (1.0 - p)
                    }
                })
                .product::<f64>()
        })
        .collect();
    let bursty: Vec<f64> = noise
        .iter()
        .zip(&weight)
        .map(|(z, w)| z * w.sqrt())
        .collect();

    let q = default_q();
    let scales = default_scales(n, 1);
    for (name, x) in [("plain noise", &noise), ("cascade volatility", &bursty)] {
        let spec = mfdfa(x, &q, &scales, 1).expect("long enough");
        println!(
            "{name:<19} Hurst h(2) = {:.2}   spectrum width = {:.2}",
            spec.hurst(),
            spec.width()
        );
    }
}
