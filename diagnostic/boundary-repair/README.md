# Two-sided boundary repair

This directory contains the selected experimental estimator for continuous
densities on a known finite interval, its paired validation against the
current probit estimator, and the retained endpoint-classifier ablation.

## Selected estimator

The reference implementation assumes the support has been mapped affinely to
$[0,1]$. First fit the current unweighted probit estimator $\widehat f_B$ and
use its CDF to define sample-shrinking endpoint weights. With
$r_n=\min(0.25,n^{-1/2})$ and $s(t)=3t^2-2t^3$,

$$
w_L(x)=1-s\left(\min\left\{1,\frac{\widehat F_B(x)}{r_n}\right\}\right),
\qquad
w_U(x)=1-s\left(\min\left\{1,\frac{1-\widehat F_B(x)}{r_n}\right\}\right),
$$

and $w_B(x)=1-w_L(x)-w_U(x)$. The finite-boundary expert averages Gaussian
equivalent-kernel local-linear and local-quadratic kernels using one degree-2
bandwidth selected from all endpoint distances. Density contributions remain
unweighted, so the estimand does not change. The native evaluator computes the
expert only where its fusion weight is nonzero, truncates negative values, and
normalizes once after fusion.

Finite endpoint behavior is identified from
$k=\min(n-1,\lceil2\sqrt n\rceil)$ order statistics. At the lower endpoint,

$$
\widehat\beta_L
=
\frac{k}{\sum_{i=1}^k
\log\{X_{(k+1)}/X_{(i)}\}},
$$

with the analogous calculation on $1-X$ at the upper endpoint. The finite
expert is used if
$\widehat\beta\geq0.9$ and
$\widehat\beta(1-1.645/\sqrt{k})\leq1$; otherwise the endpoint reuses the bulk
estimator. Thus exploding, vanishing, and ambiguous endpoints require no
additional fit.

The finite or bulk endpoint components are fused with $w_L,w_B,w_U$,
truncated at zero, and normalized numerically. There is no power transform,
offset transform, oversmoothed tail expert, or cross-validation selector.

## Retained files

- `scripts/bounded_estimator.R`: pre-native tuning implementation.
- `scripts/validate_bounded_estimator.R`: paired validation against current.
- `scripts/plot_bounded_estimator_validation.R`: numerical and realization
  report.
- `scripts/ablate_endpoint_classifier.R`: paired classifier ablation.
- `scripts/plot_endpoint_classifier_ablation.R`: classifier report.
- `results/` and `plots/`: corresponding 100-replication outputs.
- `RESULTS.md`: final numerical conclusions.
- `HANDOFF.md`: concise continuation state.
- `NATIVE_CPP.md`: pointer to the shared one- and two-sided native roadmap.

Run from the R-package root:

```sh
Rscript inst/include/kde1d-cpp/diagnostic/boundary-repair/scripts/validate_bounded_estimator.R
Rscript inst/include/kde1d-cpp/diagnostic/boundary-repair/scripts/plot_bounded_estimator_validation.R
Rscript inst/include/kde1d-cpp/diagnostic/boundary-repair/scripts/ablate_endpoint_classifier.R 100
Rscript inst/include/kde1d-cpp/diagnostic/boundary-repair/scripts/plot_endpoint_classifier_ablation.R
```
