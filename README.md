# mfdfa

[![crates.io](https://img.shields.io/crates/v/mfdfa.svg)](https://crates.io/crates/mfdfa)
[![docs.rs](https://docs.rs/mfdfa/badge.svg)](https://docs.rs/mfdfa)

Multifractal detrended fluctuation analysis for Rust. It measures how the scaling of a time series changes with the size of its moves. You get the full MFDFA spectrum (`h(q)`, `tau(q)`, `f(alpha)`), the ordinary DFA Hurst exponent, the Higuchi fractal dimension and a wavelet Hurst estimate. There are no required dependencies. Results are checked against independent Python packages.

## Install

```toml
[dependencies]
mfdfa = "0.1"
```

## Example

```rust
use mfdfa::{default_q, default_scales, higuchi_fd, mfdfa, wavelet_hurst};

fn main() -> Result<(), mfdfa::Error> {
    // Toy returns from a seeded generator. Use your own data here.
    let mut s = 42u64;
    let returns: Vec<f64> = (0..4096)
        .map(|_| {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((s >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * 0.02
        })
        .collect();

    let spec = mfdfa(&returns, &default_q(), &default_scales(returns.len(), 1), 1)?;
    println!("DFA Hurst h(2) = {:.3}", spec.hurst());   // about 0.5 for noise
    println!("spectrum width = {:.3}", spec.width());   // about 0 for monofractal noise

    // Higuchi works on the path, not the increments.
    let mut p = 0.0;
    let path: Vec<f64> = returns.iter().map(|r| { p += r; p }).collect();
    println!("Higuchi FD     = {:.3}", higuchi_fd(&path, 10)); // about 1.5
    println!("wavelet Hurst  = {:.3}", wavelet_hurst(&returns)?);
    Ok(())
}
```

More in `examples/`: `spectrum` and `streaming` (the `ta` feature).

## What is in it

| Function | Answers | Feed it |
|---|---|---|
| `mfdfa(x, q, scales, order)` | generalised Hurst `h(q)`, `tau(q)`, singularity spectrum `f(alpha)` | increments (returns) |
| `dfa_hurst(x, scales, order)` | ordinary DFA Hurst exponent, `h(2)` | increments |
| `higuchi_fd(x, k_max)` | Higuchi fractal dimension of a curve | the path (prices) |
| `wavelet_hurst(x)` | Abry-Veitch Hurst exponent on the db4 wavelet | increments |
| `hurst::hurst_rs(prices)` | rescaled-range Hurst exponent | prices |
| `RollingHurst`, `RollingHiguchi` | the same, one bar at a time | prices or increments |

Errors are values, not panics. Short input, bad scales or NaN input return an `Error` that says what went wrong.

Features:

- `ta`: `ta::Next`, `ta::Reset` and `ta::Period` for the rolling estimators.
- `serde`: `Serialize` and `Deserialize` on every result type.
- `synth`: exact fractional Gaussian noise (`synth::fgn`, `synth::fbm`).
- `cli`: the `mfdfa` command line tool, for example `mfdfa prices.csv --column close --prices`.

Also included: Petrosian, Katz and box-counting dimensions, rescaled-range Hurst, and HP, Baxter-King, SSA and EMD detrending.

## Docs

- API docs: [docs.rs/mfdfa](https://docs.rs/mfdfa)
- Crate guide with the test story: [docs/CRATE.md](docs/CRATE.md)
- Formulas: [docs/ALGORITHMS.md](docs/ALGORITHMS.md)
- Where the code came from: [docs/PROVENANCE.md](docs/PROVENANCE.md)

## License

MIT. Carved out of [srfm-lab](https://github.com/Mattbusel/srfm-lab) (MIT). See [LICENSE](LICENSE) and [NOTICE](NOTICE).
