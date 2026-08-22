# Boundary diagnostics

The diagnostic sweep evaluates deterministic samples from distributions with
known densities and CDFs. It covers finite and one-sided supports, support
scales from `1e-4` to `1e4`, polynomial degrees 0--2, grid sizes 100--1,000,
and fixed or automatically selected bandwidths. Evaluation points approach the
boundary on a logarithmic scale.

Generate the raw CSV data and an optional summary and PDF with:

```sh
diagnostic/run_boundary_diagnostics.sh > boundary-diagnostics.csv
Rscript diagnostic/plot_boundary_diagnostics.R boundary-diagnostics.csv
```

The executable is deliberately not registered with CTest. Its output is for
comparing numerical behavior between revisions; stable invariants discovered
through the sweep should be promoted to focused regression tests.
