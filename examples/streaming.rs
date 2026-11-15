//! Rolling Hurst and Higuchi as `ta` indicators.
//!
//! cargo run --example streaming --features ta

use mfdfa::{RollingHiguchi, RollingHurst};
use ta::Next;

fn main() {
    let mut hurst = RollingHurst::new(128).unwrap();
    let mut fd = RollingHiguchi::new(128, 8).unwrap();
    let mut price = 100.0;
    let mut s = 7u64;
    for bar in 0..400 {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let r = ((s >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * 0.02;
        // A trending stretch in the middle.
        price *= 1.0
            + r
            + if (150..250).contains(&bar) {
                0.003
            } else {
                0.0
            };
        let h = hurst.next(price);
        let d = fd.next(price);
        if bar % 50 == 0 && !h.is_nan() {
            println!("bar {bar:>3}  price {price:>8.2}  {hurst} = {h:.3}  {fd} = {d:.3}");
        }
    }
}
