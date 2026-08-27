## Boundary-repair estimator

This PR changes continuous fits with known finite support in two steps.

For a one-sided support, write $d(x)$ for distance from the finite endpoint,
let $s$ be the weighted median observed distance, and set
$\epsilon=10^{-5}s$. The transformed bulk estimator uses an endpoint-anchored,
scaled member of the Box-Cox family with power parameter $\lambda=1/4$:

$$
T\{d(x)\}
=4[\{\frac{d(x)+\epsilon}{s}\}^{1/4}
-\{\frac{\epsilon}{s}\}^{1/4}],
\qquad
|T'(x)|
=\frac{1}{s}\{\frac{d(x)+\epsilon}{s}\}^{-3/4}.
$$

This Box-Cox transform replaces the previous log transform and makes the
construction scale equivariant. Two-sided fits retain the regularized probit
bulk transform.

Each finite endpoint is then classified independently from its first
$k\simeq2\sqrt{n_e}$ order statistics. In the unweighted case,

$$
\widehat\beta
=\frac{k}{\sum_{i=1}^k
\log\{d_{(k+1)}/d_{(i)}\}}.
$$

Endpoints consistent with a finite nonzero limiting density receive a boundary
expert. Exploding, vanishing, ambiguous, and degenerate endpoints retain the
bulk. The expert is a Gaussian local-linear equivalent-kernel estimate with an
automatically selected bandwidth on the original scale. The selector uses all
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

> **Before merging:** remove `diagnostic/`; it is retained in this branch only
> to make the PR validation reproducible during review.

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
| one-sided | 0 | 1.810 (0.116) | 1.134 (0.049) | 0.586 (0.015) | 0.519 (0.011) |
| one-sided | 1 | 0.789 (0.038) | 0.508 (0.011) | 0.566 (0.042) | 0.612 (0.037) |
| one-sided | 2 | 0.550 (0.034) | 0.394 (0.015) | 0.231 (0.008) | 0.231 (0.008) |
| two-sided | 0 | 0.946 (0.005) | 0.940 (0.005) | 0.943 (0.006) | 0.926 (0.006) |
| two-sided | 1 | 0.856 (0.007) | 0.817 (0.007) | 0.846 (0.005) | 0.850 (0.004) |
| two-sided | 2 | 0.881 (0.013) | 0.906 (0.006) | 0.943 (0.007) | 0.932 (0.007) |

### Boundary ratio, averaged over families

| support | degree | $n=25$ | $n=100$ | $n=1000$ | $n=2000$ |
|---|---:|---:|---:|---:|---:|
| one-sided | 0 | 3.149 (0.278) | 1.553 (0.097) | 0.563 (0.020) | 0.461 (0.014) |
| one-sided | 1 | 0.931 (0.073) | 0.437 (0.017) | 0.519 (0.047) | 0.590 (0.042) |
| one-sided | 2 | 0.612 (0.048) | 0.396 (0.024) | 0.190 (0.009) | 0.190 (0.008) |
| two-sided | 0 | 0.918 (0.005) | 0.939 (0.005) | 0.955 (0.006) | 0.940 (0.006) |
| two-sided | 1 | 0.817 (0.009) | 0.821 (0.007) | 0.878 (0.004) | 0.876 (0.004) |
| two-sided | 2 | 0.877 (0.012) | 0.910 (0.005) | 0.958 (0.006) | 0.946 (0.006) |

Against the unrepaired new bulk, averaged global ISE ratios range from
0.83--0.94 for degree 0, 0.73--0.88 for degree 1, and 0.87--0.92 for degree 2
on one-sided support. The corresponding bounded ratios are 0.93--0.95,
0.82--0.86, and 0.88--0.94.

<details>
<summary>Family-level global / boundary ISE ratios</summary>

#### One-sided, degree 0

| family | $n=25$ | $n=100$ | $n=1000$ | $n=2000$ |
|---|---:|---:|---:|---:|
| Exponential | 1.48 / 2.67 | 1.36 / 2.41 | 0.58 / 0.56 | 0.51 / 0.45 |
| Gamma$(0.75)$ | 2.05 / 2.76 | 0.91 / 0.97 | 0.33 / 0.32 | 0.28 / 0.27 |
| Gamma$(2)$ | 0.79 / 0.66 | 0.84 / 0.75 | 0.87 / 0.80 | 0.86 / 0.77 |
| half-normal | 0.88 / 1.28 | 0.79 / 1.08 | 0.54 / 0.53 | 0.46 / 0.38 |
| Lomax$(2)$ | 3.84 / 8.37 | 1.77 / 2.56 | 0.61 / 0.60 | 0.50 / 0.44 |

#### One-sided, degree 1

| family | $n=25$ | $n=100$ | $n=1000$ | $n=2000$ |
|---|---:|---:|---:|---:|
| Exponential | 0.70 / 0.83 | 0.47 / 0.43 | 0.37 / 0.25 | 0.35 / 0.24 |
| Gamma$(0.75)$ | 0.72 / 0.88 | 0.33 / 0.30 | 0.32 / 0.30 | 0.36 / 0.34 |
| Gamma$(2)$ | 0.74 / 0.59 | 0.84 / 0.73 | 0.93 / 1.03 | 1.00 / 1.21 |
| half-normal | 0.61 / 0.65 | 0.50 / 0.40 | 0.47 / 0.29 | 0.44 / 0.26 |
| Lomax$(2)$ | 1.17 / 1.71 | 0.40 / 0.32 | 0.74 / 0.72 | 0.91 / 0.90 |

#### One-sided, degree 2

| family | $n=25$ | $n=100$ | $n=1000$ | $n=2000$ |
|---|---:|---:|---:|---:|
| Exponential | 0.43 / 0.53 | 0.32 / 0.38 | 0.10 / 0.07 | 0.10 / 0.07 |
| Gamma$(0.75)$ | 0.47 / 0.56 | 0.18 / 0.17 | 0.12 / 0.11 | 0.10 / 0.10 |
| Gamma$(2)$ | 0.61 / 0.48 | 0.84 / 0.83 | 0.67 / 0.59 | 0.67 / 0.61 |
| half-normal | 0.37 / 0.42 | 0.32 / 0.32 | 0.16 / 0.10 | 0.13 / 0.06 |
| Lomax$(2)$ | 0.87 / 1.09 | 0.30 / 0.28 | 0.11 / 0.08 | 0.15 / 0.11 |

#### Two-sided, degree 0

| family | $n=25$ | $n=100$ | $n=1000$ | $n=2000$ |
|---|---:|---:|---:|---:|
| Beta$(0.75,0.75)$ | 1.01 / 1.02 | 1.03 / 1.04 | 1.05 / 1.05 | 1.06 / 1.06 |
| Beta$(0.75,2)$ | 0.98 / 0.97 | 1.01 / 1.01 | 1.09 / 1.08 | 1.07 / 1.07 |
| Beta$(1,2)$ | 0.90 / 0.89 | 0.84 / 0.85 | 0.82 / 0.85 | 0.77 / 0.80 |
| Beta$(2,0.75)$ | 0.99 / 0.98 | 1.01 / 1.01 | 1.11 / 1.10 | 1.07 / 1.06 |
| Beta$(2,1)$ | 0.94 / 0.92 | 0.87 / 0.88 | 0.76 / 0.80 | 0.78 / 0.81 |
| Beta$(2,2)$ | 0.86 / 0.68 | 0.95 / 0.88 | 1.00 / 1.00 | 1.00 / 1.00 |
| Uniform | 0.94 / 0.96 | 0.86 / 0.91 | 0.76 / 0.81 | 0.72 / 0.78 |

#### Two-sided, degree 1

| family | $n=25$ | $n=100$ | $n=1000$ | $n=2000$ |
|---|---:|---:|---:|---:|
| Beta$(0.75,0.75)$ | 0.90 / 0.90 | 0.86 / 0.88 | 0.95 / 0.95 | 0.98 / 0.98 |
| Beta$(0.75,2)$ | 0.86 / 0.86 | 0.90 / 0.90 | 0.97 / 0.97 | 0.99 / 0.99 |
| Beta$(1,2)$ | 0.79 / 0.74 | 0.71 / 0.73 | 0.72 / 0.75 | 0.69 / 0.72 |
| Beta$(2,0.75)$ | 0.87 / 0.87 | 0.88 / 0.89 | 0.96 / 0.96 | 0.99 / 0.99 |
| Beta$(2,1)$ | 0.84 / 0.81 | 0.72 / 0.75 | 0.69 / 0.74 | 0.71 / 0.73 |
| Beta$(2,2)$ | 0.93 / 0.78 | 0.97 / 0.90 | 1.00 / 1.00 | 1.00 / 1.00 |
| Uniform | 0.79 / 0.76 | 0.68 / 0.70 | 0.62 / 0.76 | 0.61 / 0.73 |

#### Two-sided, degree 2

| family | $n=25$ | $n=100$ | $n=1000$ | $n=2000$ |
|---|---:|---:|---:|---:|
| Beta$(0.75,0.75)$ | 0.90 / 0.91 | 0.94 / 0.95 | 1.04 / 1.04 | 1.07 / 1.07 |
| Beta$(0.75,2)$ | 0.87 / 0.88 | 0.95 / 0.96 | 1.09 / 1.08 | 1.08 / 1.07 |
| Beta$(1,2)$ | 0.74 / 0.77 | 0.82 / 0.83 | 0.81 / 0.83 | 0.77 / 0.79 |
| Beta$(2,0.75)$ | 0.89 / 0.89 | 0.95 / 0.96 | 1.11 / 1.10 | 1.07 / 1.07 |
| Beta$(2,1)$ | 0.89 / 0.89 | 0.80 / 0.82 | 0.76 / 0.80 | 0.77 / 0.81 |
| Beta$(2,2)$ | 1.04 / 0.96 | 1.04 / 0.96 | 1.00 / 1.00 | 1.00 / 1.00 |
| Uniform | 0.85 / 0.85 | 0.84 / 0.88 | 0.80 / 0.85 | 0.77 / 0.81 |

</details>

## Bias, variance, classification, and scale

The repair primarily reduces variance, with deliberate boundary bias. Against
pre-PR at scale 1, aggregate global variance ratios are 0.38--0.62 one-sided
for degree 0, 0.39--0.56 for degree 1, and 0.20--0.24 for degree 2. Bounded
ratios are 0.96--0.98, 0.84--0.87, and 0.87--0.96, respectively. This matches
the visual goal: substantially calmer one-sided endpoints and moderate
smoothing on bounded support.

Classifier behavior becomes appropriately conservative with sample size. For
one-sided finite targets, selection rates at $n=100/2000$ are 58%/78%
(Exponential), 75%/82% (half-normal), and 56%/76% (Lomax). False finite
selection falls from 22% to 4% for Gamma$(0.75)$ and from 38% to 0% for
Gamma$(2)$. For Uniform, either endpoint is selected in 90%/98% and both in
51%/71% of samples. Beta$(2,2)$ falls from 4.5% both-endpoint selection at
$n=100$ to 0% at $n\geq1000$.

For Exponential samples at $n=100$, the 5%--95% range of the estimated value
at the endpoint is 0.64--1.34 when the expert is selected, compared with
0.01--49.6 when the classifier retains the raw bulk. Thus the expert delivers
the intended visual stabilization when selected; remaining extreme endpoint
values are primarily classifier false negatives.


<!-- Attach boundary-examples.png here. The plotted realizations are selected
     by median candidate ISE within their intended classifier branch. -->

## Coarse runtime

Median end-to-end R fit times for the default degree-2 estimator are below.

| case | $n$ | pre-PR | candidate | slowdown |
|---|---:|---:|---:|---:|
| one active endpoint | 100 | 0.48 ms | 0.79 ms | 1.66x |
| one active endpoint | 1,000 | 1.52 ms | 1.90 ms | 1.25x |
| one active endpoint | 10,000 | 12.47 ms | 13.23 ms | 1.06x |
| two active endpoints | 100 | 0.51 ms | 0.83 ms | 1.63x |
| two active endpoints | 1,000 | 1.46 ms | 1.77 ms | 1.22x |
| two active endpoints | 10,000 | 10.97 ms | 11.40 ms | 1.04x |
| one-sided fallback | 100 / 1,000 / 10,000 | -- | -- | 1.08x / 1.03x / 1.00x |
| two-sided fallback | 100 / 1,000 / 10,000 | -- | -- | 1.00x / 0.98x / 0.99x |

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
