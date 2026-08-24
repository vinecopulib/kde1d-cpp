# Validation results

## Final estimator

The final two-boundary estimator combines the current probit bulk fit with an
equal local-linear/local-quadratic expert at endpoints classified as finite.
All other endpoint types reuse bulk. The validation uses 100 paired
replications for Uniform, Beta$(1,2)$, Beta$(0.75,2)$, and
Beta$(0.75,0.75)$ at $n=100,1000,2000$. The primary measure is integrated
squared error (ISE).

Global ISE ratios to the current bounded estimator are:

| density | $n=100$ | $n=1000$ | $n=2000$ |
|---|---:|---:|---:|
| Uniform | 0.884 | 0.830 | 0.814 |
| Beta$(1,2)$ | 0.879 | 0.818 | 0.802 |
| Beta$(0.75,2)$ | 0.943 | 1.050 | 1.045 |
| Beta$(0.75,0.75)$ | 0.945 | 1.037 | 1.033 |

Mean ratios across the four densities are 0.913, 0.934, and 0.923. The
finite-endpoint cases improve clearly. Losses for exploding densities are
small and accepted in exchange for keeping a generic, shape-agnostic finite
expert.

All retained estimates are nonnegative and normalized numerically. The
maximum mass error in the validation is $2.2\times10^{-16}$.

## Endpoint classifier

The selected rule uses one tail-index estimate with
$k=\min(n-1,\lceil2\sqrt n\rceil)$ and a 95% finite-endpoint decision. It was
compared with single scales $\sqrt n$ and $1.5\sqrt n$, and multiscale median
and voting rules. Mean ISE ratios to current are:

| $n$ | $\sqrt n$, 95% | $2\sqrt n$, 95% | multiscale vote, 90% |
|---:|---:|---:|---:|
| 100 | 0.885 | 0.912 | 0.919 |
| 1000 | 0.976 | 0.936 | 0.996 |
| 2000 | 0.994 | 0.921 | 0.981 |

The smaller scale is slightly better at $n=100$, but $2\sqrt n$ is clearly
better at $n=1000,2000$ and substantially reduces the worst exploding-tail
loss. Multiscale classification adds complexity without improving ISE or the
representative fits.

## Bulk decision for nonfinite endpoints

A final focused comparison used asymmetric Beta$(a,2)$ and symmetric
Beta$(a,a)$ densities for $a=0.9,0.75,0.6$. An oversmoothed probit tail reduced
integrated variance by roughly one third, but its mean ISE ratios to simply
reusing bulk were 0.948, 1.079, and 1.095 at $n=100,1000,2000$. Bulk won all
six large-sample scenarios. The additional transformed fit, offset, bandwidth
multiplier, and fusion branch were therefore removed from the selected
estimator.

The retained CSVs contain global and endpoint-local ISE, squared bias,
integrated variance, selection rates, minimum density, and mass error. The
PDFs add representative full-density and logarithmic-distance endpoint plots.
