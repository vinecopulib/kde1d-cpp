# Native boundary-expert implementation roadmap

## Scope

Implement both selected finite-boundary repairs through one concise internal
path:

| support | unchanged bulk estimator | candidate endpoints | bandwidth floor |
|---|---|---|---|
| exactly one finite endpoint | power-$3/4$ transform | the finite endpoint | closest 75% of boundary distances |
| two finite endpoints | probit transform | lower and upper, classified independently | all boundary distances |

The existing transformed estimator remains the bulk component in both cases.
Unbounded continuous, discrete, and zero-inflated behavior remains unchanged.

The production estimator activates for continuous fits of degrees zero through
two with at least 16 observations and effective sample size at least 16. Public
case weights and manual bandwidths are supported. Apply the bandwidth
multiplier to every expert bandwidth after selection and flooring.

Expose boundary repair as a constructor option that defaults to enabled. When
disabled, retain the transformed bulk estimator. Do not add another fitted-
object representation.

## Minimal integration point

Call the repair coordinator only when `type_ == VarType::continuous`, at least
one support bound is finite, and at least 16 cleaned observations remain. The
coordinator additionally retains the bulk fit when the effective sample size is
below 16.

Keep the current body of `Kde1d::fit()` through construction of the normalized
bulk `InterpolationGrid`. Preserve the cleaned observations on their original
scale before applying `boundary_transform()`. Immediately after constructing
the bulk grid:

1. classify every finite endpoint from original-scale boundary distances;
2. return to the unchanged bulk path if no endpoint is confidently finite;
3. fit a boundary component only for each selected endpoint, using the existing
   original-scale grid;
4. form endpoint weights from the normalized bulk CDF;
5. fuse bulk and selected endpoint values and replace `grid_` with the
   normalized fused grid.

The existing PDF, CDF, quantile, simulation, support truncation, and
log-likelihood code then operates on the fused grid. `bandwidth_` continues to
report the bulk bandwidth.

Boundary distances and local densities can stay in the original units. Every
classifier ratio is dimensionless, and the boundary bandwidths
scale with the data, so this is equivalent to the R reference's affine map to
$[0,1]$ without adding a second normalization path.

No changes should be required in `boundary_transform()`, `boundary_correct()`,
`construct_grid_points()`, `InterpolationGrid`, `PluginBandwidthSelector`, or
the public constructors.

## Shared private helpers in `kde1d.hpp`

Keep the implementation in the existing header and add three focused helpers.

### 1. Endpoint classifier

Given ascending distances from one endpoint, use
$k=\min(n-1,\lceil2\sqrt n\rceil)$ and

$$
\widehat\beta=
\frac{k}{\sum_{i=1}^k
\log\{d_{(k+1)}/\max(d_{(i)},\epsilon d_{(k+1)})\}}}.
$$

The scale-relative numerical floor preserves scale equivariance. Select the
finite expert when

$$
\widehat\beta\geq0.9,
\qquad
\widehat\beta\left\{1-
\frac{\Phi^{-1}(0.95)}{\sqrt{k}}\right\}\leq1.
$$

Run this helper once for the one-sided finite endpoint and independently on
lower and upper distances for two-sided support. Exploding, vanishing, and
ambiguous endpoints all select bulk. A nonpositive $(k+1)$st distance or a
nonfinite/nonpositive logarithmic denominator is a degenerate classifier input
and must also select bulk.

### 2. Local-polynomial boundary component

The helper takes boundary distances, ascending evaluation distances, rank
weights, a floor fraction, and degree one or two.

- Select the weighted bandwidth with `PluginBandwidthSelector`.
- Select one unweighted degree-2 bandwidth from the closest
  $\max(4,\lceil q n\rceil)$ distances, where $q=0.75$ for one-sided support
  and $q=1$ for two-sided support.
- Share it between the degree-one and degree-two kernels, then apply the public
  multiplier. Evaluate density contributions with equal observation weights
  so the estimand remains the original density.

Construct the equivalent-kernel coefficients as the first row of the inverse
truncated-normal moment matrix. Generate the moments with

$$
\mu_0(a)=\Phi(a),\qquad
\mu_1(a)=-\phi(a),\qquad
\mu_j(a)=(j-1)\mu_{j-2}(a)-a^{j-1}\phi(a).
$$

Use the explicit degree-one and degree-two determinants to avoid a dynamic
matrix inverse at every evaluation point. Average the two kernels before
fusion, clamp negative values, and return both the density and its influence
numerator for the EDF calculation described below.

Evaluate directly only where the endpoint weight is nonzero. With ordered
distances, omit observations more than six bandwidths from an evaluation
point; the omitted Gaussian tail mass is below $2\times10^{-9}$. Normalize the
fused grid once at the end.

### 3. Boundary-repair coordinator

This helper owns endpoint distances, classification, local component calls,
CDF weights, reflection, and fusion.

For an upper endpoint, compute on ascending reflected distances and reverse the
component once. The floor fraction is the only local-component policy that
differs between one- and two-sided support.

The helper should return an empty result when every endpoint selects bulk.
`Kde1d::fit()` can then leave the already constructed bulk grid untouched.

## Weights and fusion

Let $F_B(x)$ be the bulk CDF and

$$
q_n=\min(0.25,n^{-1/2}),
\qquad
s(r)=3r^2-2r^3.
$$

For a lower endpoint use

$$
w_L(x)=1-s\left(\min\{1,F_B(x)/q_n\}\right),
$$

and for an upper endpoint use

$$
w_U(x)=1-s\left(\min\{1,(1-F_B(x))/q_n\}\right).
$$

Set the weight of every bulk-classified endpoint to zero. The final values are

$$
\widehat f(x)=
w_L(x)\widehat f_L(x)
+\{1-w_L(x)-w_U(x)\}\widehat f_B(x)
+w_U(x)\widehat f_U(x).
$$

With $q_n\leq0.25$, the endpoint weights do not overlap. The same expression
covers one-sided support by setting the absent endpoint weight to zero. Clamp
the fused values at zero and construct a normalized `InterpolationGrid` on the
existing original-scale grid.

## Conditional EDF approximation

The current EDF is already a conditional influence-trace approximation: it
holds the selected bandwidth and final grid normalization fixed. Extend the
same convention to expert selection and fusion rather than returning `NaN`.

For a degree-$p$ boundary component, let $M_p(a)$ be the truncated-normal
moment matrix above, let $h$ be the shared bandwidth, and let
$\widetilde f_p$ denote its value before clamping. Its influence numerator is

$$
\nu_p(x)=
\mathbf 1\{\widetilde f_p(x)>0\}
\frac{\phi(0)}{n h}
\left[M_p\{d(x)/h\}^{-1}\right]_{00}.
$$

Only the constant equivalent-kernel coefficient remains because a
self-contribution has kernel argument zero. For the equal linear/quadratic
expert, use

$$
f_E(x)=\frac{f_1(x)+f_2(x)}{2},
\qquad
\nu_E(x)=\frac{\nu_1(x)+\nu_2(x)}{2}.
$$

For the bulk component, the existing influence column is its log-density
leverage $\ell_B(x)$, so its influence numerator is
$\nu_B(x)=f_B(x)\ell_B(x)$. Fuse the numerators with exactly the same weights
as the densities:

$$
\begin{aligned}
g(x)&=w_L(x)f_L(x)
 +\{1-w_L(x)-w_U(x)\}f_B(x)+w_U(x)f_U(x),\\
\nu(x)&=w_L(x)\nu_L(x)
 +\{1-w_L(x)-w_U(x)\}\nu_B(x)+w_U(x)\nu_U(x),\\
\ell(x)&=\frac{\nu(x)}{g(x)}.
\end{aligned}
$$

Thus the mixture leverage is weighted by local density shares, not just by
the fusion weights. Estimate

$$
\operatorname{edf}=\sum_{i=1}^n
\min\{3,\max(0,\ell(X_i))\},
$$

retaining the current per-observation cap. The final common normalization of
$g$ and $\nu$ cancels in their ratio when its derivative is ignored.

This approximation conditions on bandwidths, endpoint classifications, CDF
weights, clamp states, and normalization constants. It deliberately omits the
derivatives of those data-adaptive choices. That keeps it consistent with the
existing EDF convention and costs only additional grid-vector arithmetic; no
refits or per-observation matrices are needed. Validate it diagnostically
against the corresponding leave-one-out approximation

$$
1+\sum_{i=1}^n
\{\log \widehat f(X_i)-\log \widehat f_{-i}(X_i)\}.
$$

First confirm that this diagnostic reproduces the existing bulk EDF to the
expected approximation error. A weight-derivative correction is worth
considering only if the active-expert comparison then shows an additional
systematic discrepancy.

## Weights and remaining fit metadata

Expert density contributions use the normalized public case weights, matching
the weighted bulk estimand. Effective sample size determines the minimum sample
check, fusion width, tail order, and classifier confidence. The weighted
classifier uses a conservative lower threshold of 0.975. The one-sided
bandwidth uses the closest 75% of positive-weight observations by row count and
passes their case weights to the plug-in selector. The influence numerator uses
the local average case weight on the FFT grid.

Recompute `loglik_` from the final grid using the existing code. Use the
conditional mixture trace above for `edf_` whenever a finite expert is active,
and retain the current calculation unchanged for pure-bulk fits.

No classifier choice, tail index, boundary bandwidth, or weight state needs to
persist in the fitted object for the first implementation. Diagnostic values
belong in tests and validation scripts rather than the public class layout.

## Exact production changes

The intended implementation footprint is:

1. `include/kde1d/kde1d.hpp`
   - preserve original-scale cleaned observations;
   - add the three private helpers;
   - call the coordinator after bulk-grid construction;
   - carry the bulk and expert influence numerators through fusion for EDF.
2. `test/test.cpp`
   - deterministic golden results and fallback behavior.
3. `test/test_numerical_invariants.cpp`
   - active-expert normalization, reflection, scale, support, and quantile
     invariants.
4. `diagnostic/one-sided-repair/` and `diagnostic/boundary-repair/`
   - thin R evaluators and final baseline-versus-native reports.

Avoid modifying `interpolation.hpp`, `dpik.hpp`, `stats.hpp`, constructors, or
public getters unless implementation evidence shows that it is necessary.

## Implementation sequence

1. Synchronize this branch with the isolated PR #34 fix so
   `boundary_offset_` is `NaN` before fitting. Record the current production
   head and run the native R evaluators once to create one- and two-sided
   baseline outputs before changing estimator behavior.
2. Add the eligibility check and preserve original-scale observations. With
   the coordinator returning bulk, verify tight numerical parity for every
   existing support type.
3. Implement the shared endpoint classifier and deterministic finite/bulk
   tests for lower and upper endpoints.
4. Implement the shared degree-one/two equivalent-kernel evaluator. Validate
   both bandwidth-floor policies against the retained R references.
5. Activate one-sided fusion with the 75% floor and verify the selected
   one-sided scenarios.
6. Activate two-sided fusion with the full-sample floor and verify independent
   lower/upper decisions, including asymmetric densities.
7. Complete log-likelihood and conditional-EDF handling, numerical regression
   tests, and fit-time benchmarks.
8. Run the larger R validation against the recorded native baselines, retain
   only final comparison artifacts, and prepare the PR report.

Each step should be a separately testable commit. Avoid refactoring unrelated
fit logic while implementing the method.

## Weighted and manual-bandwidth resolution

Manual bandwidths control only the transformed bulk fit. Boundary bandwidths
remain independently selected in original data units, exactly as when the bulk
bandwidth is automatic. The public multiplier remains the common smoothing
control and scales both bulk and expert bandwidths.

Weighted experiments rejected defining the one-sided 75% bandwidth subset by
cumulative weight: concentrated weights made that subset too narrow and
unstable. Counting positive-weight observations was both simpler and more
accurate. Effective sample size gave a small but consistent improvement for the
fusion width under concentrated weights. A weighted classifier threshold of
0.975 avoided most false repair of exploding targets while retaining the
finite-boundary gains. Global weight rescaling and insertion or removal of
zero-weight observations are numerically invariant.

Integer-weight replication equivalence is not required here because the
existing bulk bandwidth selector interprets public weights through effective
sample size rather than as frequency counts. Supporting frequency weights
would require a separate change to the bulk estimator's weight semantics.

### Degrees zero and one

The same averaged local-linear/local-quadratic boundary expert is used for all
three bulk degrees. It outperformed matching-order and local-linear-only
alternatives in global and boundary ISE. The gains persisted for weighted and
manual-bandwidth fits, and boundary scale equivariance is tested for every
degree. Degree zero primarily trades substantially lower boundary bias for a
small one-sided variance increase; degree one reduces both bias and variance.

## Native tests

Add compact tests to the existing executables.

Shared coverage:

- nonnegativity, unit mass, monotone CDF, quantile round trips, finite
  log-likelihood and EDF, and support truncation;
- exact EDF parity with the current path when no expert is active, plus
  diagnostic comparison with leave-one-out log-likelihood optimism when an
  expert is active;
- scaling over several orders of magnitude;
- ties, observations exactly at an endpoint, small samples, and nearly
  degenerate samples;
- regression coverage for unbounded, discrete, zero-inflated, weighted,
  fixed-bandwidth, and nondefault-degree fits;
- fit-time benchmarks at representative $n$ values.

One-sided coverage:

- deterministic Exponential, half-normal, Gamma$(2)$, Gamma$(0.75)$, and
  Lomax samples;
- lower/upper reflection while the expert is active;
- finite and bulk classifier outcomes;
- golden grid or PDF values against `one_sided_estimator.R`.

Two-sided coverage:

- deterministic Uniform, Beta$(1,2)$, Beta$(0.75,2)$, and
  Beta$(0.75,0.75)$ samples;
- lower and upper finite/bulk outcomes, including asymmetric decisions;
- affine equivariance on non-unit intervals;
- golden grid or PDF values against `bounded_estimator.R`.

If direct local evaluation is materially expensive, optimize only after a
profile identifies it as the relevant fit cost.

## R validation without C++ experiment bloat

Keep Monte Carlo comparisons in R. Use small native-evaluator scripts that run
the installed estimator for deterministic scenario seeds and write
replication-level global and boundary ISE plus aggregate squared bias and
variance. Run them once at the recorded baseline commit and once at the
candidate commit, then join the outputs in reporting scripts.

This avoids a legacy-mode switch in the public C++ API and avoids maintaining
two production estimator paths solely for benchmarking.

The final one-sided suite should expand the retained scenarios with
intermediate sample sizes, steeper finite-endpoint slopes, heavier remote
tails, upper-bounded reflections, and visual realizations. The final two-sided
suite should retain varied symmetric and asymmetric finite, vanishing, and
exploding endpoint shapes plus local boundary plots. Report paired ISE ratios
with Monte Carlo uncertainty, global and local bias/variance, classifier rates
computed by matching R code, visual examples, and fit time relative to the
current transformed baselines.
