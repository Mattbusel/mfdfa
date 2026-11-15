# Algorithms

Exact definitions of what each function computes, so results can be compared
with other implementations.

## MFDFA (`mfdfa`, `fluctuation_function`, `dfa_hurst`)

Kantelhardt, Zschiegner, Koscielny-Bunde, Havlin, Bunde and Stanley,
"Multifractal detrended fluctuation analysis of nonstationary time series",
Physica A 316 (2002) 87-114.

Input: a series `x_1..x_N` of increments (for prices, log returns).

1. Profile: `Y(i) = sum_{k=1..i} (x_k - mean(x))`.
2. For each scale `s`, `N_s = floor(N / s)` segments from the start
   (`Y[v s .. (v+1) s]`) and `N_s` from the end (`Y[N mod s + v s ..]`), so
   `2 N_s` in total.
3. In each segment, fit a least-squares polynomial of degree `order` against
   the sample index and keep the mean squared residual `F^2(v, s)`. The
   polynomial basis is orthonormalised (Gram-Schmidt on a scaled abscissa), so
   high orders stay stable.
4. `F_q(s) = ( (1 / 2N_s) sum_v [F^2(v, s)]^{q/2} )^{1/q}` for `q != 0`, and
   `F_0(s) = exp( (1 / 4N_s) sum_v ln F^2(v, s) )`.
5. `h(q)` is the OLS slope of `ln F_q(s)` on `ln s`.
6. `tau(q) = q h(q) - 1`, `alpha(q) = d tau / d q`, `f(alpha) = q alpha - tau(q)`.
   The derivative is `numpy.gradient(tau, q)`: second-order central
   differences inside (exact for non-uniform `q` spacing), first-order at the
   two ends.

`dfa_hurst` is step 5 for `q = 2` alone.

Constraints: `q` strictly increasing (at least two values), scales distinct,
each at least `order + 2`, and `N >= 2 max(scales)`.

`default_q()` is -5 to 5 step 0.5; `default_scales(N, order)` is up to 20
log-spaced integers from `max(10, 2 (order + 1))` to `N / 4`.

Step 4 differs from the `MFDFA` Python package only in that the package drops
`q = 0`; on every other `q` the two agree to 1e-10 (`tests/reference.rs`).

## Higuchi fractal dimension (`higuchi`, `higuchi_fd`)

Higuchi, "Approach to an irregular time series on the basis of the fractal
theory", Physica D 31 (1988) 277-283.

Input: the curve `x_1..x_N` itself.

For `k = 1..k_max` and `m = 1..k`, with `M = floor((N - m) / k)`:

```text
L_m(k) = (N - 1) / (k^2 M) * sum_{j=1..M} |x_{m + j k} - x_{m + (j-1) k}|
L(k)   = mean over m of L_m(k)
D      = OLS slope of ln L(k) on ln(1 / k)
```

Needs `k_max >= 2` and `N >= 2 k_max`. `higuchi_fd` returns NaN where
`higuchi` returns an error. The estimate is not clamped. `antropy.higuchi_fd`
adds `1e-9` to its regression denominator; with that removed the two agree to
1e-10.

## Wavelet Hurst (`wavelet_hurst`, `wavelet_hurst_detail`)

Abry and Veitch, "Wavelet analysis of long-range-dependent traffic", IEEE
Transactions on Information Theory 44 (1998) 2-15.

The orthonormal Daubechies db4 filter (8 taps, PyWavelets `dec_lo`/`dec_hi`)
with periodic boundaries: coefficient `i` at a level is
`sum_k h[k] x[(2 i + 4 - k) mod n]`, which is PyWavelets'
`mode="periodization"` alignment. If a level has odd length its last sample
is dropped before the next split.

For octaves `j` with at least 8 coefficients, `y_j = log2( mean(d_j^2) )`; the
slope of `y_j` on `j` is fitted by least squares with weights `n_j` (number of
coefficients), and `H = (slope + 1) / 2`. This is for increments: for fGn of
exponent `H`, `E[d_j^2]` grows like `2^{j (2H - 1)}`.

## Rescaled range (`hurst::hurst_rs`)

Kept from srfm-lab: log returns of the prices, R/S averaged over
non-overlapping blocks of 8, 16, 32, ... up to `N / 2` plus the full length,
slope of `ln(R/S)` on `ln(n)`, clamped to [0, 1]. No Anis-Lloyd small-sample
correction, so short windows read high.

## Fractional Gaussian noise (`synth::fgn`)

Davies and Harte, "Tests for Hurst effect", Biometrika 74 (1987) 95-101.
The autocovariance `g(k) = (|k+1|^{2H} - 2|k|^{2H} + |k-1|^{2H}) / 2` is
embedded in a circulant of size `2N`; its eigenvalues come from one FFT, and
the real part of the FFT of `sqrt(lambda_k / 2N) (Z_k + i Z'_k)` is an exact
sample. The Gaussian draws come from SplitMix64 and Box-Muller, so a seed
reproduces the series on every platform.
