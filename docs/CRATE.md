# mfdfa

Multifractal detrended fluctuation analysis (MFDFA), Higuchi fractal dimension
and wavelet Hurst estimation for time series, in Rust, with no required
dependencies.

```toml
[dependencies]
mfdfa = "0.1"
```

| Function | Answers | Feed it |
|---|---|---|
| `mfdfa(x, q, scales, order)` | generalised Hurst `h(q)`, `tau(q)`, singularity spectrum `f(alpha)` | increments (returns) |
| `dfa_hurst(x, scales, order)` | ordinary DFA Hurst exponent, `h(2)` | increments |
| `higuchi_fd(x, k_max)` | Higuchi fractal dimension of a curve | the path (prices, log prices) |
| `wavelet_hurst(x)` | Abry-Veitch Hurst exponent on the db4 wavelet | increments |
| `hurst::hurst_rs(prices)` | rescaled-range Hurst exponent | prices |
| `RollingHurst`, `RollingHiguchi` | the same, one bar at a time | prices or increments |

## MFDFA in five lines

```rust
use mfdfa::{mfdfa, default_q, default_scales};

let returns: Vec<f64> = load_returns();          // your data
let spec = mfdfa(&returns, &default_q(), &default_scales(returns.len(), 1), 1)?;

println!("DFA Hurst h(2) = {:.3}", spec.hurst());
println!("spectrum width = {:.3}", spec.width());  // ~0 for monofractal noise
for (q, h) in spec.q.iter().zip(&spec.h_q) {
    println!("h({q:+.1}) = {h:.3}");
}
```

`MfdfaSpectrum` holds `q`, `scales`, the full fluctuation function `fq[q][s]`,
`h_q`, `tau_q = q h(q) - 1`, `alpha` and `f_alpha` (one entry per `q`; the
Legendre derivative uses the same finite differences as `numpy.gradient`).

Errors are values, not panics: too-short input, scales smaller than
`order + 2`, NaN input, or a constant series with `q <= 0` all return an
`Error` saying which.

## Higuchi and wavelet Hurst

```rust
let d = mfdfa::higuchi_fd(&log_prices, 10);   // ~1.5 for a random walk, NaN if undefined
let h = mfdfa::wavelet_hurst(&returns)?;       // ~0.5 for uncorrelated returns
```

For fractional Brownian motion `D` is about `2 - H`, so the two estimates
cross-check each other.

## Streaming, with the `ta` crate

```toml
mfdfa = { version = "0.1", features = ["ta"] }
```

```rust
use mfdfa::{RollingHurst, RollingHiguchi};
use ta::Next;

let mut hurst = RollingHurst::new(128)?;        // fed prices, DFA-1 on log returns
let mut fd = RollingHiguchi::new(128, 8)?;      // fed prices
for bar in bars {
    let h = hurst.next(&bar);                   // any T: ta::Close, or a plain f64
    let d = fd.next(bar.close);
    // NaN until the window is full
}
```

Both also implement `ta::Reset` and `ta::Period`, and work without the feature
through `update(value)`.

## Features

| Feature | Adds |
|---|---|
| `ta` | `ta::Next<f64>`, `ta::Next<&T: Close>`, `ta::Reset`, `ta::Period` for the rolling estimators |
| `serde` | `Serialize` / `Deserialize` on every result type |
| `synth` | `synth::fgn(n, hurst, seed)` and `synth::fbm`: exact fractional Gaussian noise (Davies-Harte, uses `rustfft`) |
| `cli` | the `mfdfa` binary: `mfdfa prices.csv --column close --prices` |

## How it is checked

`tests/reference.rs` rebuilds two test signals bit for bit and compares against
values produced by independent Python packages (`tools/gen_reference.py`):

- `F_q(s)`, `h(q)`, `tau(q)`, `alpha`, `f(alpha)` against the
  [MFDFA](https://github.com/LRydin/MFDFA) package, detrending orders 1 to 3,
  to 1e-10;
- Higuchi dimension against [antropy](https://github.com/raphaelvallat/antropy)
  for five `k_max` values, to 1e-10;
- db4 coefficients against PyWavelets `mode="periodization"`, to 1e-12.

`tests/known_answers.rs` checks closed forms: the binomial p-model cascade's
`h(q)` (Kantelhardt et al. 2002), and DFA, wavelet and Higuchi estimates on
exact fractional Gaussian noise at H = 0.3, 0.5, 0.7, 0.85.

More: [ALGORITHMS.md](https://github.com/Mattbusel/mfdfa/blob/main/docs/ALGORITHMS.md) for the definitions,
[PROVENANCE.md](https://github.com/Mattbusel/mfdfa/blob/main/docs/PROVENANCE.md) for what changed from srfm-lab.

## Also in the box

Ported from srfm-lab with their tests: Petrosian, Katz and box-counting
dimensions (`fractal_dim_ext`, `fractal_dimension`), rescaled-range Hurst
(`hurst`), HP, Baxter-King, SSA and EMD detrending (`detrending`), and the
price-oriented `multifractal::MultifractalAnalyzer`.

## License

MIT. Carved out of [srfm-lab](https://github.com/Mattbusel/srfm-lab) (MIT); see
NOTICE.
