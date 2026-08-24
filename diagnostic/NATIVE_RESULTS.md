# Native boundary-expert evidence

The final comparison with the pre-PR estimator has not yet been run. Its
protocol is in [`README.md`](README.md). In particular, the one-sided primary
benchmark will be the actual production log-transform fit from pre-PR `main`,
not the unrepaired fourth-root fit used during development.

The following implementation checks remain relevant to the selected method.

## Selected implementation

- Finite endpoints use the equal average of local-linear and local-quadratic
  Gaussian equivalent kernels; the same expert is used for bulk degrees zero,
  one, and two.
- Both kernels share an automatically selected degree-two bandwidth. Selection
  uses all positive-weight observations on bounded support and the closest 75%
  by row count on one-sided support. A manual bandwidth controls only the bulk;
  the public multiplier scales both bandwidths.
- Public case weights are used by the classifier, expert estimand, expert
  bandwidth selector, effective sample size, and influence approximation.
- A cubic CDF-based fusion region has probability width
  $\min(0.25,n_e^{-1/2})$. Endpoints classified as exploding, vanishing, or
  ambiguous retain the bulk fit.

## FFT accuracy and runtime

The endpoint kernels are pointwise combinations of an ordinary Gaussian KDE
and its first two derivatives, evaluated on a shared 256-bin FFT grid. In the
100-replication development suite, the mean relative ISE perturbation from FFT
binning was at most 0.011% in every support/sample-size group. Among active
experts, the median absolute paired ISE change was 0.023%, the 99th percentile
was 0.18%, and the maximum was 0.66%.

The direct-kernel fallback was rejected. Across 51--401 output grid points at
$n=25,100$, always using FFT slowed the worst case by 1.25x. On the default
grid at $n=100$, one active endpoint changed from 0.64 to 0.67 ms and two
active endpoints from 0.77 to 0.74 ms. At $n=2000$, FFT reduced one-endpoint
fit time from 1.21 to 0.92 ms and two-endpoint time from 2.25 to 0.96 ms.

## EDF approximation

The conditional mixture EDF was checked against

$$
1+\sum_i\{\log\widehat f(X_i)-\log\widehat f_{-i}(X_i)\}.
$$

For deterministic finite-endpoint samples at $n=25,50,100$, discrepancies
were 0.08--0.30 EDF units for Uniform and Exponential examples. For
Beta$(1,2)$ they were 0.27--0.33 at $n=50,100$ and 0.84 at $n=25$. This was
considered adequate for an asymptotic trace approximation that already holds
bandwidths fixed. Pure-bulk fits retain the prior EDF calculation exactly.
