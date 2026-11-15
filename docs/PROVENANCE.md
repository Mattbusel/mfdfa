# Provenance: what came from srfm-lab and what changed

mfdfa started as the `fractal-engine` crate in
[srfm-lab](https://github.com/Mattbusel/srfm-lab) (MIT), commit
`8cbfcb0551d44e3ed1c20bdae9dc805d61521c62` (2026-11-07), with
`lib/math/fractal.py` as the Python cross-check. This file lists what was
kept, what was fixed and what was left behind, so anyone comparing the two
codebases knows why numbers differ.

## New

- `mfdfa()`: general MFDFA with any detrending order (srfm-lab's Rust version
  was linear only and took prices), `q = 0`, and full `alpha` / `f(alpha)`
  arrays. Checked against the `MFDFA` Python package.
- `higuchi()` / `higuchi_fd()` at the crate root: unclamped, any length of at
  least `2 k_max`. Checked against `antropy`.
- `wavelet_hurst()`: srfm-lab had a DWT but no wavelet Hurst estimator.
- `RollingHurst`, `RollingHiguchi` with `ta` trait implementations.
- `synth::fgn` / `synth::fbm` for known-answer tests.

## Fixed

| Where (srfm-lab) | Problem | Now |
|---|---|---|
| `lib/math/fractal.py::higuchi_fd` | Leaves out the `1 / k` factor in `L_m(k)`, so the slope is `D - 1`, then clips to [1, 2]. It returns about 1.0 for a random walk (antropy: 1.50) and about 1.0 for white noise (antropy: 2.00). | Full Higuchi normalisation, matches antropy. (The Rust `fractal_dimension.rs` version had the factor and was right.) |
| `wavelet.rs::inverse_dwt` | Not an inverse: a 64-sample signal came back with a maximum error of 2.6. The unit test only checked the output length. | Exact transpose of the orthonormal periodic transform; reconstruction to 1e-10 is tested for lengths 64, 100 and 257. |
| `wavelet.rs` filters | db4 taps at 15 digits, time-reversed relative to PyWavelets, with a different phase. | Full-precision PyWavelets taps and alignment; coefficients match `pywt.dwt(..., mode="periodization")`. |
| `multifractal.rs` | `h(q)` clamped to [0, 2] and zero-variance segments silently dropped, which changes `F_q(s)` for negative `q`. Spectrum built from interior points only. | No clamping; a zero-variance segment with `q <= 0` is reported as `Error::Degenerate`. Spectrum has one point per `q`. |
| `hurst.rs` test `rs_single_constant_returns_none` | Asserted `... || true`, so it could not fail. | Asserts the constant case returns `None`; an oscillating-series test now uses the previously dead helper. |

## Kept with their tests

`hurst.rs` (rescaled range, `HurstRegime`), `fractal_dimension.rs` (box
counting, `FdClass`; the Higuchi wrapper is renamed `higuchi_analysis` and
calls the new core), `fractal_dim_ext.rs` (Petrosian, Katz, dispatcher,
`FdTradingRegime`), `detrending.rs` (HP, Baxter-King, SSA, EMD) and
`multifractal.rs` (`MultifractalAnalyzer` is now a wrapper over `mfdfa()`).
The 30 srfm-lab integration tests run unchanged apart from the crate name, as
`tests/legacy_modules.rs`. Strategy-specific names (LARSA, "Black Hole") were
removed from doc comments.

## Left behind

- `main.rs`: the srfm-lab CLI with a hard-coded `data/crypto_trades.csv`
  default; replaced by the `cli` feature binary.
- `pattern_library.rs`, `similarity.rs`, `regime_detector.rs`: trading
  heuristics built on top of the estimators, not estimators themselves.
- Dependencies `ndarray`, `statrs`, `anyhow` (unused by the kept modules) and
  the always-on `serde`, `serde_json`, `clap`, `csv`. `rustfft` is kept, behind
  the `synth` feature, for exact fGn synthesis.
