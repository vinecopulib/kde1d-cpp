## Boundary-repair estimator

This PR changes continuous fits with known finite support in two steps.

For a one-sided support, write $d(x)$ for distance from the finite endpoint,
let $s$ be the weighted median observed distance, and set
$\epsilon=10^{-5}s$. The transformed bulk estimator now uses

$$
T\{d(x)\}
=4\left[\left\{\frac{d(x)+\epsilon}{s}\right\}^{1/4}
-\left(\frac{\epsilon}{s}\right)^{1/4}\right],
\qquad
|T'(x)|
=\frac{1}{s}\left\{\frac{d(x)+\epsilon}{s}\right\}^{-3/4}.
$$

This replaces the previous log transform and makes the construction scale
equivariant. Two-sided fits retain the regularized probit bulk transform.

Each finite endpoint is then classified independently from its first
$k\simeq2\sqrt{n_e}$ order statistics. In the unweighted case,

$$
\widehat\beta
=\frac{k}{\sum_{i=1}^k
\log\{d_{(k+1)}/d_{(i)}\}}.
$$

Endpoints consistent with a finite nonzero limiting density receive a boundary
expert. Exploding, vanishing, ambiguous, and degenerate endpoints retain the
bulk. The expert is a Gaussian local-linear equivalent-kernel estimate,

using an automatically selected degree-2 plug-in bandwidth on the original
scale. The selector uses all
positive-weight observations on bounded support and the closest 75% by row
count on one-sided support. The two required Gaussian convolutions are
evaluated on a shared 256-bin FFT grid.

Let $F_B$ be the normalized bulk CDF,

$$
q=\min(1/4,n_e^{-1/2}),
\qquad
g(p)=\{1-\min(1,p/q)^2\}^2.
$$

Selected lower and upper experts receive the biweight fusion weights
$w_L(x)=g\{F_B(x)\}$ and $w_U(x)=g\{1-F_B(x)\}$, giving

$$
\widehat f(x)
=w_L(x)\widehat f_L(x)
+\{1-w_L(x)-w_U(x)\}\widehat f_B(x)
+w_U(x)\widehat f_U(x).
$$

The result is clamped at zero and normalized. The same boundary expert is used
for bulk degrees 0, 1, and 2. EDF is a conditional influence-trace
approximation that holds the classifier, selected bandwidths, fusion weights,
clamps, and normalization fixed.

## Validation design

The large-scale tables below predate the biweight and local-linear
simplifications and will be refreshed before merge.

The primary benchmark is the actual pre-PR package at parent commit `7faec53`
with C++ commit `b16d514`: log transform for one-sided support and probit for
two-sided support. The experiment used 200 paired replications for
$n\in\{25,100,1000,2000\}$ and degrees 0, 1, and 2.

- One-sided: Exponential, half-normal, Gamma$(0.75)$, Gamma$(2)$, and
  Lomax$(2)$, with lower/upper reflections and scales $10^{-4},1,10^4$.
- Two-sided: Uniform and Beta$(1,2)$, Beta$(2,1)$, Beta$(0.75,2)$,
  Beta$(2,0.75)$, Beta$(0.75,0.75)$, and Beta$(2,2)$.
- Global ISE and boundary ISE were integrated separately. The boundary region
  is the outer 10% probability mass at each finite endpoint.

The tables below use scale 1; reflected one-sided results agree to numerical
precision. Entries are candidate/baseline ratios, so **values below 1 favor the
PR**. Parentheses are paired Monte Carlo standard errors.

### Global ISE ratio, averaged over families

| support | degree | $n=25$ | $n=100$ | $n=1000$ | $n=2000$ |
|---|---:|---:|---:|---:|---:|
| one-sided | 0 | 0.748 (0.044) | 0.579 (0.028) | 0.513 (0.011) | 0.475 (0.010) |
| one-sided | 1 | 0.505 (0.020) | 0.411 (0.011) | 0.684 (0.063) | 0.724 (0.050) |
| one-sided | 2 | 0.211 (0.036) | 0.198 (0.021) | 0.208 (0.018) | 0.260 (0.015) |
| two-sided | 0 | 0.951 (0.005) | 0.943 (0.005) | 0.941 (0.006) | 0.925 (0.005) |
| two-sided | 1 | 0.867 (0.006) | 0.823 (0.007) | 0.847 (0.005) | 0.852 (0.004) |
| two-sided | 2 | 0.892 (0.013) | 0.911 (0.006) | 0.940 (0.006) | 0.932 (0.006) |

### Boundary ratio, averaged over families

| support | degree | $n=25$ | $n=100$ | $n=1000$ | $n=2000$ |
|---|---:|---:|---:|---:|---:|
| one-sided | 0 | 0.702 (0.049) | 0.513 (0.030) | 0.447 (0.013) | 0.399 (0.011) |
| one-sided | 1 | 0.399 (0.022) | 0.303 (0.011) | 0.636 (0.075) | 0.702 (0.059) |
| one-sided | 2 | 0.144 (0.027) | 0.151 (0.027) | 0.159 (0.016) | 0.205 (0.015) |
| two-sided | 0 | 0.923 (0.005) | 0.944 (0.005) | 0.954 (0.005) | 0.939 (0.005) |
| two-sided | 1 | 0.824 (0.008) | 0.826 (0.007) | 0.879 (0.004) | 0.877 (0.004) |
| two-sided | 2 | 0.883 (0.012) | 0.913 (0.005) | 0.955 (0.005) | 0.945 (0.005) |

The expert also improves over the unrepaired new bulk. Averaged global ISE
ratios range from 0.84--0.94 for degree 0, 0.74--0.89 for degree 1, and
0.88--0.92 for degree 2 on one-sided support. The corresponding bounded ratios
are 0.93--0.95, 0.82--0.87, and 0.89--0.94.

<details>
<summary>Family-level global / boundary ISE ratios</summary>

#### One-sided, degree 0

| family | $n=25$ | $n=100$ | $n=1000$ | $n=2000$ |
|---|---:|---:|---:|---:|
| Exponential | 0.75 / 0.73 | 0.49 / 0.45 | 0.37 / 0.30 | 0.34 / 0.27 |
| Gamma$(0.75)$ | 0.55 / 0.54 | 0.51 / 0.50 | 0.45 / 0.44 | 0.44 / 0.43 |
| Gamma$(2)$ | 0.71 / 0.52 | 0.76 / 0.58 | 0.81 / 0.69 | 0.80 / 0.67 |
| half-normal | 0.46 / 0.38 | 0.41 / 0.32 | 0.31 / 0.21 | 0.28 / 0.17 |
| Lomax$(2)$ | 1.28 / 1.34 | 0.72 / 0.70 | 0.63 / 0.59 | 0.51 / 0.45 |

#### One-sided, degree 1

| family | $n=25$ | $n=100$ | $n=1000$ | $n=2000$ |
|---|---:|---:|---:|---:|
| Exponential | 0.49 / 0.37 | 0.31 / 0.22 | 0.38 / 0.26 | 0.38 / 0.27 |
| Gamma$(0.75)$ | 0.27 / 0.25 | 0.25 / 0.21 | 0.45 / 0.43 | 0.53 / 0.52 |
| Gamma$(2)$ | 0.69 / 0.46 | 0.81 / 0.63 | 0.93 / 0.98 | 1.01 / 1.19 |
| half-normal | 0.45 / 0.32 | 0.38 / 0.23 | 0.43 / 0.25 | 0.43 / 0.25 |
| Lomax$(2)$ | 0.62 / 0.60 | 0.31 / 0.22 | 1.23 / 1.26 | 1.27 / 1.28 |

#### One-sided, degree 2

| family | $n=25$ | $n=100$ | $n=1000$ | $n=2000$ |
|---|---:|---:|---:|---:|
| Exponential | 0.07 / 0.04 | 0.02 / 0.01 | 0.06 / 0.04 | 0.09 / 0.06 |
| Gamma$(0.75)$ | 0.03 / 0.03 | 0.04 / 0.03 | 0.23 / 0.21 | 0.23 / 0.21 |
| Gamma$(2)$ | 0.38 / 0.18 | 0.72 / 0.56 | 0.59 / 0.44 | 0.62 / 0.50 |
| half-normal | 0.09 / 0.04 | 0.03 / 0.02 | 0.03 / 0.01 | 0.06 / 0.03 |
| Lomax$(2)$ | 0.49 / 0.43 | 0.17 / 0.13 | 0.13 / 0.09 | 0.30 / 0.23 |

#### Two-sided, degree 0

| family | $n=25$ | $n=100$ | $n=1000$ | $n=2000$ |
|---|---:|---:|---:|---:|
| Beta$(0.75,0.75)$ | 1.02 / 1.03 | 1.03 / 1.04 | 1.04 / 1.04 | 1.05 / 1.05 |
| Beta$(0.75,2)$ | 0.99 / 0.98 | 1.01 / 1.00 | 1.07 / 1.07 | 1.06 / 1.06 |
| Beta$(1,2)$ | 0.93 / 0.91 | 0.85 / 0.86 | 0.83 / 0.86 | 0.78 / 0.81 |
| Beta$(2,0.75)$ | 1.00 / 0.99 | 1.01 / 1.00 | 1.09 / 1.08 | 1.06 / 1.05 |
| Beta$(2,1)$ | 0.95 / 0.93 | 0.89 / 0.89 | 0.78 / 0.81 | 0.79 / 0.82 |
| Beta$(2,2)$ | 0.82 / 0.65 | 0.94 / 0.88 | 1.00 / 1.00 | 1.00 / 1.00 |
| Uniform | 0.95 / 0.98 | 0.88 / 0.93 | 0.77 / 0.82 | 0.74 / 0.79 |

#### Two-sided, degree 1

| family | $n=25$ | $n=100$ | $n=1000$ | $n=2000$ |
|---|---:|---:|---:|---:|
| Beta$(0.75,0.75)$ | 0.91 / 0.91 | 0.87 / 0.88 | 0.94 / 0.95 | 0.97 / 0.97 |
| Beta$(0.75,2)$ | 0.88 / 0.87 | 0.89 / 0.90 | 0.96 / 0.97 | 0.98 / 0.98 |
| Beta$(1,2)$ | 0.81 / 0.76 | 0.72 / 0.74 | 0.73 / 0.76 | 0.69 / 0.73 |
| Beta$(2,0.75)$ | 0.89 / 0.88 | 0.88 / 0.88 | 0.95 / 0.96 | 0.98 / 0.98 |
| Beta$(2,1)$ | 0.86 / 0.82 | 0.74 / 0.76 | 0.70 / 0.75 | 0.71 / 0.74 |
| Beta$(2,2)$ | 0.91 / 0.73 | 0.96 / 0.90 | 1.00 / 1.00 | 1.00 / 1.00 |
| Uniform | 0.82 / 0.78 | 0.70 / 0.72 | 0.64 / 0.77 | 0.62 / 0.74 |

#### Two-sided, degree 2

| family | $n=25$ | $n=100$ | $n=1000$ | $n=2000$ |
|---|---:|---:|---:|---:|
| Beta$(0.75,0.75)$ | 0.91 / 0.92 | 0.95 / 0.96 | 1.02 / 1.02 | 1.05 / 1.05 |
| Beta$(0.75,2)$ | 0.89 / 0.89 | 0.94 / 0.94 | 1.06 / 1.06 | 1.06 / 1.05 |
| Beta$(1,2)$ | 0.77 / 0.80 | 0.84 / 0.85 | 0.82 / 0.85 | 0.78 / 0.81 |
| Beta$(2,0.75)$ | 0.90 / 0.90 | 0.94 / 0.95 | 1.08 / 1.07 | 1.05 / 1.05 |
| Beta$(2,1)$ | 0.90 / 0.90 | 0.83 / 0.85 | 0.78 / 0.82 | 0.79 / 0.83 |
| Beta$(2,2)$ | 0.99 / 0.89 | 1.02 / 0.95 | 1.00 / 1.00 | 1.00 / 1.00 |
| Uniform | 0.87 / 0.87 | 0.86 / 0.90 | 0.82 / 0.86 | 0.79 / 0.83 |

</details>

## Bias, variance, classification, and scale

The ISE gains are primarily a variance reduction, with some deliberate
small-sample boundary bias. Against pre-PR at scale 1, aggregate global
variance ratios are 0.27--0.53 one-sided for degree 0, 0.29--0.43 for degree 1,
and 0.04--0.12 for degree 2. Bounded ratios are 0.99--1.03, 0.86--0.90, and
0.89--1.00, respectively. This matches the visual goal: substantially calmer
one-sided endpoints and moderate smoothing on bounded support.

Classifier behavior becomes appropriately conservative with sample size. For
one-sided finite targets, selection rates at $n=100/2000$ are 58%/78%
(Exponential), 75%/82% (half-normal), and 56%/76% (Lomax). False finite
selection falls from 22% to 4% for Gamma$(0.75)$ and from 38% to 0% for
Gamma$(2)$. For Uniform, either endpoint is selected in 90%/98% and both in
51%/71% of samples. Beta$(2,2)$ falls from 4.5% both-endpoint selection at
$n=100$ to 0% at $n\geq1000$.

For Exponential samples at $n=100$, the 5%--95% range of the estimated value
at the endpoint is 0.56--1.44 when the expert is selected, compared with
0.01--49.6 when the classifier retains the raw bulk. Thus the expert delivers
the intended visual stabilization when selected; remaining extreme endpoint
values are primarily classifier false negatives.


<!-- Attach boundary-examples.png here. The plotted realizations are selected
     by median candidate ISE within their intended classifier branch. -->

## Coarse runtime

Median end-to-end R fit times for the default degree-2 estimator are below.

| case | $n$ | pre-PR | candidate | slowdown |
|---|---:|---:|---:|---:|
| one active endpoint | 100 | 0.48 ms | 0.81 ms | 1.70x |
| one active endpoint | 1,000 | 1.51 ms | 1.96 ms | 1.29x |
| one active endpoint | 10,000 | 12.43 ms | 13.57 ms | 1.09x |
| two active endpoints | 100 | 0.51 ms | 0.87 ms | 1.71x |
| two active endpoints | 1,000 | 1.45 ms | 1.85 ms | 1.27x |
| two active endpoints | 10,000 | 10.97 ms | 11.83 ms | 1.08x |
| one-sided fallback | 100 / 1,000 / 10,000 | -- | -- | 1.09x / 1.05x / 1.02x |
| two-sided fallback | 100 / 1,000 / 10,000 | -- | -- | 1.00x / 1.01x / 1.01x |

## Less common paths and numerical checks

- Public case weights are used consistently by the classifier, expert
  estimand, bandwidth selector, effective sample size, and EDF approximation.
  Tests cover weight rescaling and zero-weight-row invariance.
- A manually supplied bandwidth controls the transformed bulk only. Boundary
  experts retain their original-scale automatic bandwidth; the public
  multiplier scales both. This unusual path has focused regression tests but
  was not included in the full Monte Carlo sweep.
- Setting `boundary_repair=false` recovers the transformed bulk and skips all
  expert work.
- Tests cover degrees 0--2, normalization, nonnegativity, monotone CDFs,
  quantile round trips, support truncation, reflection, scales $10^{-4}$ to
  $10^4$, weights, fixed bandwidths, and pure-bulk EDF parity.
- The conditional EDF approximation differed from leave-one-out optimism by
  about 0.08--0.30 EDF units for deterministic Uniform/Exponential examples;
  the largest checked discrepancy was 0.84 for Beta$(1,2)$ at $n=25$.
