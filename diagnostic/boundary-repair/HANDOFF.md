# Two-boundary repair handoff

## Current state

The directory now contains one experimental estimator for continuous data
with two known finite bounds. It fits the current probit bulk estimator once,
adds an equal local-linear/local-quadratic expert for confidently finite
endpoints, and otherwise reuses bulk. Smooth CDF-scale gates shrink as
$n^{-1/2}$. The finite-endpoint classifier uses
$k=\min(n-1,\lceil2\sqrt n\rceil)$ order statistics and a 95% decision.

There is no exploding-tail expert, power warp, transformation offset, or CV
selector. Superseded oversmoothing and offset code and outputs were removed.
The retained 100-replication validation and classifier ablation are summarized
in `RESULTS.md`.

## Reproduce

From the R-package root:

```sh
Rscript inst/include/kde1d-cpp/diagnostic/boundary-repair/scripts/validate_bounded_estimator.R
Rscript inst/include/kde1d-cpp/diagnostic/boundary-repair/scripts/plot_bounded_estimator_validation.R
Rscript inst/include/kde1d-cpp/diagnostic/boundary-repair/scripts/ablate_endpoint_classifier.R 100
Rscript inst/include/kde1d-cpp/diagnostic/boundary-repair/scripts/plot_endpoint_classifier_ablation.R
```

For a smoke test, pass a smaller replication count and temporary output path
as the first two arguments to either simulation script.

## Known limitations

- The R reference assumes support $[0,1]$; native code must apply an affine
  map for arbitrary finite bounds.
- Weighted observations have not been specified or validated for the tail
  index, weighted ranks, or shrinking gate.
- The finite-versus-bulk decision is hard, so between-sample endpoint behavior
  can still change discretely.
- Before release, test affine support changes, weaker and stronger endpoint
  singularities, interior modes, ties, and weighted data.

## Next work

The immediate next investigation is the analogous one-sided estimator. Keep
that work in a separate sibling diagnostic directory so this two-boundary
reference remains stable. Native implementation of this estimator should
follow `NATIVE_CPP.md` after the one-sided design has settled enough to assess
shared infrastructure.

Unrelated top-level generated documentation changes were left untouched.
