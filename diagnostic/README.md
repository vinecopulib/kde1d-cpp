# Boundary-expert validation

This directory contains the retained validation path for the native estimator.
The estimator itself, including the support transforms, endpoint classifier,
boundary kernels, fusion, bandwidth handling, and EDF approximation, is
documented in [`../docs/overview-continuous.dox`](../docs/overview-continuous.dox).

## Final comparison

`validate_native_boundary_experts.R` evaluates one installed implementation at
a time. Run it with these method labels:

| method | installed revision | estimator |
|---|---|---|
| `pre-pr` | pre-PR `main` (`1260b5e` for the planned comparison) | production log transform on one-sided support and probit transform on two-sided support |
| `bulk` | candidate | fourth-root/probit bulk with `boundary_repair = FALSE` |
| `expert` | candidate | fourth-root/probit bulk with the selected boundary experts |

Using the actual pre-PR revision is intentional: emulating its log transform in
the candidate would miss other grid and numerical behavior of the production
benchmark.

The suite covers lower- and upper-bounded versions of five one-sided densities
at scales `1e-4`, `1`, and `1e4`, seven two-sided densities, and sample sizes
25, 100, 1,000, and 2,000. It writes replication-level global and boundary ISE,
EDF, and log-likelihood, plus integrated squared bias and variance summaries.
Paired seeds are shared across methods.

From the R-package root, install each requested revision into a separate
library and run, for example:

```sh
Rscript inst/include/kde1d-cpp/diagnostic/validate_native_boundary_experts.R \
  pre-pr diagnostic-results/pre-pr.csv 100 /path/to/pre-pr-library
Rscript inst/include/kde1d-cpp/diagnostic/validate_native_boundary_experts.R \
  bulk diagnostic-results/bulk.csv 100 /path/to/candidate-library
Rscript inst/include/kde1d-cpp/diagnostic/validate_native_boundary_experts.R \
  expert diagnostic-results/expert.csv 100 /path/to/candidate-library
```

Generated CSVs and figures are not versioned. The final aggregate comparison
and representative plots should be recorded in `NATIVE_RESULTS.md` after the
validation is complete.

`NATIVE_RESULTS.md` currently retains only implementation-level evidence that
still applies to the selected estimator. Superseded R prototypes, candidate
estimators, tuning outputs, and the standalone C++ diagnostic were removed;
stable numerical properties are covered by the unit and invariant tests.
