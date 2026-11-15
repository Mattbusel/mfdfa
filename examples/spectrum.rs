//! Compute an MFDFA spectrum, Higuchi dimension and wavelet Hurst for a
//! synthetic multifractal-ish return series (volatility clustering).
//!
//! cargo run --example spectrum

use mfdfa::{default_q, default_scales, higuchi_fd, mfdfa, wavelet_hurst};

fn main() {
    // Returns with slowly varying volatility: a crude stand-in for real data.
    let mut s = 2024u64;
    let mut uniform = move || {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (s >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut vol = 0.01;
    let returns: Vec<f64> = (0..8192)
        .map(|_| {
            vol = (0.98 * vol + 0.02 * 0.01 + 0.002 * (uniform() - 0.5)).abs();
            vol * (uniform() - 0.5) * 3.46
        })
        .collect();

    let spec = mfdfa(&returns, &default_q(), &default_scales(returns.len(), 1), 1)
        .expect("series is long enough");
    println!("{:>6} {:>8} {:>8} {:>8}", "q", "h(q)", "alpha", "f(a)");
    for i in (0..spec.q.len()).step_by(2) {
        println!(
            "{:>6.1} {:>8.4} {:>8.4} {:>8.4}",
            spec.q[i], spec.h_q[i], spec.alpha[i], spec.f_alpha[i]
        );
    }
    println!("DFA Hurst h(2)   = {:.4}", spec.hurst());
    println!("spectrum width   = {:.4}", spec.width());

    let mut p = 0.0;
    let log_price: Vec<f64> = returns
        .iter()
        .map(|r| {
            p += r;
            p
        })
        .collect();
    println!("Higuchi FD       = {:.4}", higuchi_fd(&log_price, 10));
    println!("wavelet Hurst    = {:.4}", wavelet_hurst(&returns).unwrap());
}
