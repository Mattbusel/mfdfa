# mfdfa: one Hurst exponent hides the risk you care about

[![crates.io](https://img.shields.io/crates/v/mfdfa.svg)](https://crates.io/crates/mfdfa)
[![docs.rs](https://docs.rs/mfdfa/badge.svg)](https://docs.rs/mfdfa)

A single Hurst exponent tells you whether a series trends or mean-reverts on average. It does not tell you whether the big moves scale like the small ones, and in markets that is the part that hurts. Two return series can have the same Hurst exponent and completely different tails:

```sh
cargo run --example same_hurst
```
```text
plain noise         Hurst h(2) = 0.49   spectrum width = 0.10
cascade volatility  Hurst h(2) = 0.47   spectrum width = 0.74
```

Multifractal detrended fluctuation analysis (Kantelhardt et al. 2002) measures exactly that difference: a generalised Hurst exponent `h(q)` for every moment `q`, and the singularity spectrum `f(alpha)` whose width is near 0 for a monofractal series and grows with intermittency. If your regime or risk model stops at one Hurst number, it is reporting the average and dropping the tails.

**Why this instead of the alternatives.** In Python, MFDFA, Higuchi and wavelet Hurst usually come from three different packages. mfdfa gives you all of them in one crate with no required dependencies, and the tests compare it with the independent Python packages (`MFDFA` 0.4.3 for `F_q(s)`, `antropy` 0.2.3 for Higuchi, PyWavelets for the db4 filters), not with its own output. Errors are values: short input, bad scales or NaN return an `Error` that says which, never a panic and never a quiet NaN spectrum.

**Why today.** `cargo add mfdfa`, no build scripts, no BLAS. The `ta` feature gives rolling Hurst and Higuchi as `ta::Next` indicators, one bar at a time, and the `cli` feature reads a CSV of prices.

**Who it is for:** quants and researchers who classify regimes or size risk from scaling behaviour, and anyone porting a Python fractal-analysis notebook to Rust who needs the same numbers.

**Trade-offs, on purpose:** `f64` slices in, plain structs out; no plotting; polynomial detrending of order you choose, nothing adaptive behind your back.

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
