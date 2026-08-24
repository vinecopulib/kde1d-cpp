# One-sided boundary repair

This workspace applies the selected two-boundary idea to continuous densities
with exactly one known finite endpoint. It is deliberately separate from
`boundary-repair/`, which is stable.

## Current candidate

The power-$3/4$ transformation estimator is the unweighted bulk fit. A single
order-statistic tail index classifies the finite endpoint using
$k=\min(n-1,\lceil2\sqrt n\rceil)$ and the same 95% rule as in the
two-boundary estimator. A confidently finite endpoint receives the averaged
local-linear/local-quadratic kernel expert; exploding, vanishing, and ambiguous
endpoints reuse bulk.

The boundary expert is blended with bulk through a smooth weight based on the
bulk CDF, with boundary probability $\min(0.25,n^{-1/2})$. Both kernels share
the degree-2 bandwidth fitted to the closest 75% of boundary distances. The
native evaluator only computes expert values where this weight is nonzero and
normalizes once after fusion. This preserves the stabilizing information in
most of the sample without allowing the most remote observations to dominate.
The construction is reflection- and scale-equivariant.

The validation suite uses 30 paired replications at $n=100,1000$ for
Exponential, half-normal, Gamma$(2)$, Gamma$(0.75)$, and Lomax$(2)$ densities.
The log transform was rejected and its experimental implementation and outputs
have been removed.

## Files

- `scripts/one_sided_estimator.R`: pre-native tuning implementation.
- `scripts/ablate_local_floor_fraction.R`: selected floor-fraction ablation.
- `scripts/plot_local_floor_fraction_ablation.R`: tuning curves and example
  estimates.
- `results/local-floor-fraction-ablation.csv` and
  `plots/local-floor-fraction-ablation.pdf`: current tuning evidence.
- `RESULTS.md`: numerical conclusions.
- `HANDOFF.md`: implementation state and next steps.
- `NATIVE_CPP.md`: shared one- and two-sided native implementation roadmap.

Run from the R-package root:

```sh
Rscript inst/include/kde1d-cpp/diagnostic/one-sided-repair/scripts/ablate_local_floor_fraction.R 30
Rscript inst/include/kde1d-cpp/diagnostic/one-sided-repair/scripts/plot_local_floor_fraction_ablation.R
```
