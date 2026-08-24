# Native boundary-expert validation

## Statistical validation

The native candidate was compared with commit `4f3d148`, immediately before
the expert implementation. The paired experiment used 100 replications at
$n=25,100,1000,2000$. One-sided scenarios were Exponential, half-normal,
Gamma$(2)$, Gamma$(0.75)$, and Lomax$(2)$; two-sided scenarios were Uniform and
six symmetric or asymmetric beta densities. Ratios below average the
scenario-specific ratios of mean ISE. Parentheses contain paired Monte Carlo
standard errors.

| support | $n$ | global ISE ratio | boundary ISE ratio |
|---|---:|---:|---:|
| one-sided | 25 | 0.932 (0.008) | 0.896 (0.016) |
| one-sided | 100 | 0.929 (0.009) | 0.903 (0.016) |
| one-sided | 1000 | 0.907 (0.012) | 0.857 (0.017) |
| one-sided | 2000 | 0.881 (0.015) | 0.824 (0.020) |
| two-sided | 25 | 0.874 (0.021) | 0.862 (0.020) |
| two-sided | 100 | 0.911 (0.008) | 0.917 (0.007) |
| two-sided | 1000 | 0.939 (0.007) | 0.954 (0.006) |
| two-sided | 2000 | 0.944 (0.007) | 0.955 (0.006) |

Average integrated variance ratios range from 0.93 to 0.94 one-sided and from
0.87 to 0.96 two-sided. Boundary variance ratios range from 0.88 to 0.90 and
from 0.86 to 0.97, respectively. The extra small-sample bias is the intended
price for smoother boundary realizations; at $n=1000,2000$, average integrated
squared bias is also lower than baseline.

The main accepted losses remain the difficult shapes identified by the R
experiments. At $n=2000$, Gamma$(0.75)$ has global/boundary ratios
1.085/1.098, Beta$(0.75,2)$ has 1.075/1.072, and the reflected
Beta$(2,0.75)$ has 1.027/1.026. Gamma$(2)$ and Beta$(2,2)$ select bulk in every
large-sample replication and are exactly unchanged. Uniform, Beta$(1,2)$, and
Beta$(2,1)$ have $n=2000$ global ratios 0.844, 0.817, and 0.830.

The deterministic native tests cover lower/upper reflection, independent
two-sided endpoint decisions, exploding/vanishing fallback, and affine support
changes. Their golden values were refreshed for the simplified evaluator.

## Runtime simplification ablation

The direct native evaluator was simplified after the main validation. The
local-linear and local-quadratic equivalent kernels now share the degree-2
bandwidth and are averaged in one observation loop. On bounded support this
bandwidth uses all endpoint distances; on one-sided support it uses the closest
75%. Endpoint kernels are evaluated only where their smooth fusion weights are
nonzero, Gaussian terms beyond six bandwidths are omitted, and the fused grid
is normalized once at the end.

A paired 100-replication check compared this implementation with commit
`cbde530` over Exponential, half-normal, Gamma$(2)$, Gamma$(0.75)$, Lomax$(2)$,
Uniform, Beta$(1,2)$, Beta$(0.75,2)$, Beta$(0.75,0.75)$, and Beta$(2,2)$.
Ratios compare the simplified evaluator with `cbde530`; parentheses contain
paired Monte Carlo standard errors:

| support | $n$ | global ISE ratio | boundary ISE ratio |
|---|---:|---:|---:|
| one-sided | 25 | 0.996 (0.002) | 0.999 (0.004) |
| one-sided | 100 | 0.983 (0.003) | 0.982 (0.003) |
| one-sided | 1000 | 0.993 (0.003) | 0.992 (0.003) |
| one-sided | 2000 | 0.992 (0.003) | 0.991 (0.003) |
| two-sided | 25 | 1.010 (0.003) | 1.014 (0.006) |
| two-sided | 100 | 0.999 (0.002) | 1.002 (0.004) |
| two-sided | 1000 | 0.999 (0.003) | 0.992 (0.004) |
| two-sided | 2000 | 0.994 (0.003) | 0.986 (0.004) |

The largest scenario-level increases were about 3.5% for Gamma$(0.75)$ and
3.9% for Beta$(0.75,0.75)$; Uniform, Beta$(1,2)$, Exponential, half-normal,
and Lomax generally improved. This supports the simpler shared-bandwidth
kernel average without a material aggregate accuracy cost.

## EDF check

The conditional mixture EDF was compared with

$$
1+\sum_i\{\log\widehat f(X_i)-\log\widehat f_{-i}(X_i)\}.
$$

For deterministic finite-endpoint samples at $n=25,50,100$, the discrepancy
was 0.08--0.30 EDF units for Uniform and Exponential examples. The asymmetric
Beta$(1,2)$ discrepancy was 0.27--0.33 at $n=50,100$ and 0.84 at $n=25$.
This is adequate for the existing asymptotic EDF convention and does not
justify a weight-derivative correction. Pure-bulk fits retain the old EDF
exactly.

## Fit time

Release builds were timed against `4f3d148`. Bulk-classified fits incur only
sorting and remain close to baseline. With a 401-point output grid, the
simplified active-expert timings are:

| active endpoints | $n=100$ | $n=1000$ | $n=2000$ |
|---|---:|---:|---:|
| one | 0.65 ms (1.7x) | 0.94 ms (2.0x) | 1.22 ms (2.1x) |
| two | 0.78 ms (2.0x) | 1.63 ms (3.7x) | 2.26 ms (4.4x) |

Relative to `cbde530`, this is a speedup of 2.6--8.7x for one active endpoint
and 3.7--8.3x for two. Absolute expert-fit cost remains below 2.3 ms at
$n=2000$ in this benchmark.

### FFT evaluation

The endpoint kernel can be written as a pointwise combination of an ordinary
Gaussian KDE and its first two derivatives. The implementation now evaluates
these three translation-invariant convolutions on a regular 256-bin grid and
interpolates them to the endpoint grid.

Against the direct endpoint evaluator, total fit times with the default
401-point grid were:

| active endpoints | $n$ | direct | FFT | speedup |
|---|---:|---:|---:|---:|
| one | 1000 | 0.93 ms | 0.78 ms | 1.19x |
| one | 2000 | 1.21 ms | 0.92 ms | 1.32x |
| two | 1000 | 1.62 ms | 0.84 ms | 1.93x |
| two | 2000 | 2.25 ms | 0.96 ms | 2.35x |

The direct-kernel fallback was removed after a separate $n=25,100$ benchmark.
Across 51--401 grid points, the worst total-fit slowdown from always using FFT
was 1.25x. With the default grid at $n=100$, one endpoint slowed from 0.64 to
0.67 ms while two endpoints improved slightly from 0.77 to 0.74 ms. In the
100-replication validation suite, the mean relative ISE perturbation from FFT
binning was at most 0.011% for every support/sample-size group. Among fits
whose endpoint expert was active, the median absolute paired ISE change was
0.023%, the 99th percentile was 0.18%, and the maximum was 0.66%.

## Reproduction

Install the package with a clean native rebuild at each revision and run
`validate_native_boundary_experts.R`, giving each run a distinct method label
and output file. The script writes replication-level ISE/EDF/log-likelihood and
a companion bias/variance summary. Join the two raw files by support,
scenario, sample size, and replication for paired ratios and uncertainty.
