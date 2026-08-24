# One-sided boundary-repair results

## Estimator held fixed during the ablation

The candidate uses the power-$3/4$ estimator as unweighted bulk. If the
order-statistic classifier identifies a finite nonzero endpoint, equal
local-linear and local-quadratic estimates are blended into the shrinking
boundary region. Only the bandwidth-floor sample fraction changes.

The experiment uses 30 paired replications for five densities at $n=100$ and
$n=1000$. It compares the lowest 15%, 20%, 25%, 30%, 35%, 40%, 50%, 75%, and
100% of boundary distances, plus a shrinking subset of $\lceil\sqrt n\rceil$
distances. All values below are candidate/current ratios, so smaller is better.

## Aggregate result

| $n$ | floor subset | global ISE | boundary ISE | squared bias | variance |
|---:|---:|---:|---:|---:|---:|
| 100 | 25% | 0.989 | 0.986 | 0.895 | 1.007 |
| 100 | 50% | 0.933 | 0.899 | 0.860 | 0.949 |
| 100 | 75% | **0.902** | **0.854** | **0.829** | **0.917** |
| 100 | 100% | 0.941 | 0.911 | 1.004 | 0.923 |
| 100 | $\sqrt n$ | 1.073 | 1.127 | 0.905 | 1.102 |
| 1000 | 25% | 0.933 | 0.896 | 0.801 | 0.976 |
| 1000 | 50% | 0.909 | 0.863 | 0.797 | 0.946 |
| 1000 | 75% | **0.889** | **0.834** | **0.789** | **0.924** |
| 1000 | 100% | 0.983 | 0.959 | 0.963 | 0.986 |
| 1000 | $\sqrt n$ | 1.043 | 1.064 | 0.806 | 1.112 |

The old 25% floor is too local. Increasing the subset through 75% reduces both
bias and variance. A shrinking $\sqrt n$ subset is clearly too small and
reintroduces instability. Using the full sample overshoots the optimum because
remote tail observations again inflate the floor.

## Selected 75% floor by density

| density | $n$ | global ISE | boundary ISE | squared bias | variance |
|---|---:|---:|---:|---:|---:|
| Exponential | 100 | 0.808 | 0.709 | 0.575 | 0.841 |
| Exponential | 1000 | 0.805 | 0.722 | 0.455 | 0.890 |
| half-normal | 100 | 0.869 | 0.744 | 0.817 | 0.878 |
| half-normal | 1000 | 0.808 | 0.672 | 0.868 | 0.801 |
| Gamma$(2)$ | 100 | 1.026 | 1.083 | 1.174 | 1.019 |
| Gamma$(2)$ | 1000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Gamma$(0.75)$ | 100 | 0.985 | 0.983 | 0.913 | 0.992 |
| Gamma$(0.75)$ | 1000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Lomax$(2)$ | 100 | 0.819 | 0.752 | 0.670 | 0.857 |
| Lomax$(2)$ | 1000 | 0.833 | 0.774 | 0.622 | 0.926 |

The full-sample floor would improve the finite light-tail cases further, but it
raises Lomax global ISE to 1.075 at $n=100$ and 1.392 at $n=1000$. The 75%
floor retains large improvements for both light and long remote tails and is
therefore the selected default.

The realization plots show the same trade-off. Relative to 25%, the selected
floor substantially calms the near-endpoint estimates. Relative to 100%, it
avoids the flattened Lomax boundary while remaining visually smooth for
Exponential data. All simulated estimates are nonnegative and normalized to
numerical precision.
