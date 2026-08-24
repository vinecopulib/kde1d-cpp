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

The summary reports four complementary diagnostics:

- PDF and CDF error against known distributions near the boundary;
- equivariance when the same data are expressed at different scales;
- agreement between left-bounded fits and reflected right-bounded fits;
- sensitivity to the interpolation grid size.

The metrics deliberately remain separate. Relative density error is not useful
where the true density vanishes, while CDF error can remain small despite a
severely unstable pointwise density. Candidate repairs should therefore be
checked against both finite nonzero boundary densities and densities that tend
to zero.

The separate transformation-boundary investigation, its retained simulation
results, and its current recommendations live in
[`boundary-repair/`](boundary-repair/README.md).

The follow-up for supports with exactly one finite endpoint lives in
[`one-sided-repair/`](one-sided-repair/README.md).
