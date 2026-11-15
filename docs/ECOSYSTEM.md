# Where these estimators fit in other libraries

The point of carving mfdfa out of srfm-lab is that the estimators are useful
outside it. These are the places where a Hurst or Higuchi indicator is missing
today, checked against each project on 2026-11-15, and what a contribution would
look like. Each contribution stands on its own: none of them needs this crate
as a dependency except the Rust one, and each is a port of a published
definition with its own tests.

## pandas-ta-classic (Python)

- **Gap.** `pandas_ta_classic/statistics/` has `entropy`, `kurtosis`, `mad`,
  `stdev`, `zscore` and others, but no Higuchi dimension and no Hurst exponent;
  a code search for "hurst" in the repository returns nothing.
- **Their rule.** Indicators that TA-Lib and tulipy do not cover need an
  independent plain-loop reference in `tests/fixtures/reference_ports.py`,
  which `tests/oracle/test_reference_ports.py` compares against.
- **Contribution.** `statistics/higuchi.py` (rolling Higuchi FD over a
  `length` window, `kmax` parameter) plus a `higuchi_fd` reference port in
  `reference_ports.py` transcribed from Higuchi (1988). The numbers to check the
  port against are the same ones this crate's `tests/reference.rs` checks
  (antropy, epsilon removed). A rolling `hurst` (DFA-1 or R/S, Hurst 1951 /
  Peng 1994) can follow as a second PR once the first lands.

## nautilus_trader (Rust core, Python API)

- **Gap.** `crates/indicators/src` has `average`, `momentum`, `volatility`,
  `ratio` and others, but no Hurst or Higuchi indicator. The repository's own
  `hurst_vpin_directional` example strategy computes a rescaled-range Hurst
  inline (`strategy.rs`, `estimate_hurst`) because there is no indicator to
  reuse.
- **Contribution.** Open an issue first proposing a `HurstExponent` indicator
  (rolling, `Indicator` trait, PyO3 binding like the other indicators) that the
  example strategy could then use; then the PR. The indicator would be written
  in their style inside `crates/indicators`, with `RollingHurst` here as the
  tested reference for expected values.

## stockstats (Python)

- **Gap.** No match for "hurst" anywhere in the repository.
- **Contribution.** A `hurst` column (`df['hurst_100']`-style window suffix,
  like the existing indicators), DFA-1 on log returns, with a unit test against
  values from this crate or from a direct numpy transcription.

## The `ta` crate (Rust)

`RollingHurst` and `RollingHiguchi` already implement `ta::Next`, `ta::Reset`
and `ta::Period`, so users of `ta` can use them today with
`features = ["ta"]`. No upstream change needed.

## What not to bring

`fractal_dim_ext::FdTradingRegime`, `hurst::HurstRegime` and the regime
thresholds are opinionated trading rules from srfm-lab. They are fine inside
this crate, but contributions to the libraries above should ship only the
estimator and leave thresholds to users.
