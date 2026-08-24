# One-sided repair handoff

The production branch already uses a scale-equivariant power-$3/4$ transform
for one-sided supports. The diagnostic candidate keeps that estimator as the
unweighted bulk and adds the selected two-boundary method's finite expert:
an averaged local-linear/local-quadratic kernel, a shrinking CDF weight, and a
$2\sqrt n$ order-statistic finite-versus-bulk decision.

The native kernels share the degree-2 bandwidth fitted to the lowest 75% of
boundary distances. A 30-replication ablation compared fixed fractions from
15% through 100% and a $\sqrt n$ shrinking subset. Performance improved up to
75%; using all observations reversed the long-tail gains because remote Lomax
observations again controlled the floor. The 75% choice has mean global ISE
ratios to current of 0.902 at $n=100$ and 0.889 at $n=1000$, with corresponding
boundary ratios 0.854 and 0.834.

The representative plots confirm that 75% removes most residual finite-endpoint
spiking. It is visually close to the full-sample floor for Exponential and
half-normal cases, but remains stable for Lomax. The main remaining loss is a
small-sample false finite classification for Gamma$(2)$: global ISE is 1.026
and boundary ISE is 1.083 relative to current at $n=100$. At $n=1000$ the
classifier always chooses bulk for both Gamma examples, so the result is
exactly current.

The log transform was rejected and removed. The working diagnostic now contains
only the current power-transform candidate and the floor-fraction validation.

Before native C++ implementation, expand the scenario suite around the selected
75% value and classifier edge cases. In particular, add finite-boundary shapes
with steep nonzero slopes, heavier remote tails, upper-bounded reflections, and
sample sizes between 100 and 1000. Keep paired global/boundary ISE, bias,
variance, normalization, nonnegativity, reflection, scale, and realization
checks.
