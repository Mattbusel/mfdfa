# Changelog

All notable changes to this crate are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the crate uses
[semantic versioning](https://semver.org/).

## [0.1.1] - 2026-11-20

### Added
- `examples/same_hurst.rs`: two return series with the same DFA Hurst exponent
  (0.49 and 0.47) and spectrum widths of 0.10 and 0.74.

### Changed
- README opens with the stance (one Hurst exponent hides the tails), why this
  crate, why today, who it is for and the trade-offs; Cargo description matches.
  No library code change.

## [0.1.0] - 2026-11-15

First release, carved out of the `fractal-engine` crate and `lib/math/fractal.py`
in [srfm-lab](https://github.com/Mattbusel/srfm-lab) at commit `8cbfcb0`.

### Added

- `mfdfa(x, q, scales, order) -> Result<MfdfaSpectrum>`: Kantelhardt MFDFA with
  any polynomial detrending order, forward and backward segments, `q = 0` via
  the logarithmic average. `MfdfaSpectrum` carries `fq`, `h_q`, `tau_q`,
  `alpha`, `f_alpha` plus `width()`, `hurst()`, `delta_h()`, `alpha_peak()`.
- `fluctuation_function`, `dfa_hurst`, `profile`, `default_q`, `default_scales`.
- `higuchi_fd(x, k_max) -> f64` and `higuchi(x, k_max) -> Result<Higuchi>`.
- `wavelet_hurst(x)` and `wavelet_hurst_detail(x, min_octave, max_octave)`:
  Abry-Veitch estimator on the db4 DWT.
- `RollingHurst` and `RollingHiguchi` streaming estimators; with the `ta`
  feature they implement `ta::Next<f64>`, `ta::Next<&T: Close>`, `ta::Reset`
  and `ta::Period`.
- `synth` feature: exact fractional Gaussian noise and fBm (Davies-Harte).
- `serde` feature: derives on all result types.
- `cli` feature: the `mfdfa` binary (CSV column in, spectrum out).
- Reference tests against the `MFDFA`, `antropy` and `PyWavelets` Python
  packages (agreement to 1e-10), closed-form checks on the binomial cascade and
  on fractional Gaussian noise.

### Changed (relative to srfm-lab)

- The db4 DWT now uses the PyWavelets periodization alignment and
  `inverse_dwt` is the exact inverse (the srfm-lab inverse did not reconstruct
  its input).
- `MultifractalAnalyzer` is a thin wrapper over `mfdfa()`; `h(q)` is no longer
  clamped to [0, 2] and zero-variance segments are no longer dropped silently.
- The legacy Higuchi wrapper is now `fractal_dimension::higuchi_analysis`;
  `higuchi_fd` at the crate root returns the unclamped estimate.
- No required dependencies; `serde`, `clap` and `csv` are optional, `ndarray`,
  `statrs` and `anyhow` are gone.

[0.1.0]: https://github.com/Mattbusel/mfdfa/releases/tag/v0.1.0
