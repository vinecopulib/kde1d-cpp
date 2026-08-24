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
| `pre-pr` | parent `7faec53`, with C++ `b16d514` | production log transform on one-sided support and probit transform on two-sided support |
| `bulk` | candidate | fourth-root/probit bulk with `boundary_repair = FALSE` |
| `expert` | candidate | fourth-root/probit bulk with the selected boundary experts |

Using the actual pre-PR revision is intentional: emulating its log transform in
the candidate would miss other grid and numerical behavior of the production
benchmark.

The suite covers lower- and upper-bounded versions of five one-sided densities
at scales `1e-4`, `1`, and `1e4`, seven two-sided densities, and sample sizes
25, 100, 1,000, and 2,000 for bulk degrees zero, one, and two. It writes
replication-level global and boundary ISE, EDF, and log-likelihood, plus
integrated squared bias and variance summaries. Paired seeds are shared across
methods and degrees.

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
is generated with:

```sh
Rscript inst/include/kde1d-cpp/diagnostic/summarize_native_boundary_validation.R \
  diagnostic-results/pre-pr.csv diagnostic-results/bulk.csv \
  diagnostic-results/expert.csv diagnostic-results/ratios
```

`benchmark_native_boundary_experts.R` performs coarse end-to-end R timings.
`extract_native_boundary_examples.R` and
`plot_native_boundary_examples.R` recreate the representative density montage
using classifier-stratified median-ISE realizations selected by the summary
script.

Final results are recorded in `NATIVE_RESULTS.md`; `PR_COMMENT.md` contains the
copy-ready methodological description and report. Superseded R prototypes,
candidate estimators, tuning outputs, and the standalone C++ diagnostic were
removed; stable numerical properties are covered by unit and invariant tests.
