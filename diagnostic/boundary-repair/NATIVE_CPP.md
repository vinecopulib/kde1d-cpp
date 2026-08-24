# Native C++ implementation plan

## Integration point

The production estimator belongs in
`include/kde1d/kde1d.hpp`, where `Kde1d::fit()` currently handles bounded
continuous data through a single probit-transformed local-likelihood fit. The
new path should activate only when both bounds are finite and the variable is
continuous. One-sided, unbounded, discrete, and zero-inflated behavior should
remain unchanged.

Avoid implementing the method through recursively constructed public
`Kde1d` objects: a nested bounded fit would activate the boundary selector
again. Extract the existing bounded-probit calculation into a private helper,
or introduce an internal mode that requests the legacy transformation fit
directly.

## Components to implement

1. Normalize observations and grid locations to $[0,1]$ using the supplied
   finite support. Convert final density values back with the inverse support
   length. This must preserve affine equivariance.
2. Fit the unweighted legacy probit bulk estimator once. Retain its PDF, CDF,
   and final original-scale grid.
3. Compute $k=\min(n-1,\lceil2\sqrt n\rceil)$ and the two order-statistic tail
   indices. Use a partial selection algorithm if convenient, although a full
   sort is acceptable initially. Protect zero endpoint distances with the C++
   equivalent of `std::numeric_limits<double>::epsilon()`.
4. Classify each endpoint as finite or bulk using the fixed finite-endpoint
   rule in `README.md`. Exploding, vanishing, and ambiguous cases all use bulk.
5. Construct rank-based lower and upper bandwidth weights using average ranks
   for ties and boundary fraction $\min(0.25,n^{-1/2})$.
6. Add a private Gaussian equivalent-kernel evaluator for boundary polynomial
   degrees one and two. Select each bandwidth through
   `bandwidth::PluginBandwidthSelector` with endpoint weights, impose the
   corresponding full-sample unweighted bandwidth as its floor, but evaluate
   density contributions with equal observation weights. Reflect data and
   grid locations for the upper endpoint.
7. Average the degree-one and degree-two endpoint densities. Clamp negative
   component values to zero before component normalization, matching the R
   reference.
8. Obtain fusion probabilities from the bulk CDF and evaluate the smoothstep
   gates. Substitute bulk or finite component values according to each
   endpoint decision, fuse arithmetically, clamp at zero, and normalize on the
   final grid.
9. Store the fused values in `interp::InterpolationGrid` so existing PDF,
    CDF, quantile, and simulation methods continue to work without a new
    public fitted-object representation.

## API and metadata decisions

The first implementation can keep the method automatic for two finite bounds,
but an internal feature switch is useful while validating. Decide whether the
selector threshold and gate constants remain private or become advanced API
parameters; the current recommendation is to keep them fixed until broader
validation justifies exposing them.

`loglik_` can be recomputed from the final interpolated density. Effective
degrees of freedom do not combine trivially across selected and fused experts;
do not silently reuse the bulk value. Either derive the fused influence trace
or return `NaN` for `edf_` until a principled definition is implemented.

The public API accepts case weights, while the retained experiment is
unweighted. Before enabling the new estimator for weighted fits, define:

- weighted ranks and tie handling for endpoint bandwidth weights;
- an effective sample size for $k$ and the shrinking gate;
- a weighted order-statistic tail index;
- whether selector confidence adjustments use raw or effective sample size.

Until those choices are validated, weighted bounded fits should remain on the
legacy path rather than pretending the unweighted selector applies.

## Tests required

- Golden PDF vectors against `scripts/bounded_estimator.R` for fixed Uniform,
  Beta$(1,2)$, Beta$(0.75,2)$, and Beta$(0.75,0.75)$ samples.
- Explicit samples that select bulk and finite behavior at each endpoint,
  including asymmetric lower/upper decisions.
- Nonnegativity, unit mass, monotone CDF, inverse-CDF consistency, and support
  truncation.
- Affine-equivariance tests on intervals such as $[-3,7]$.
- Stability with ties, observations exactly at a bound, small samples, and
  nearly degenerate samples.
- Regression tests proving that one-sided, unbounded, discrete,
  zero-inflated, and weighted fits retain their existing behavior.
- A native rerun of the retained 100-replication validation before making the
  new path the default.
