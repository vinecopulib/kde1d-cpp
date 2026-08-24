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

The deterministic native tests agree with the retained R reference values to
within $2\times10^{-4}$ relative tolerance and cover lower/upper reflection,
independent two-sided endpoint decisions, exploding/vanishing fallback, and
affine support changes.

## EDF check

The conditional mixture EDF was compared with

$$
1+\sum_i\{\log\widehat f(X_i)-\log\widehat f_{-i}(X_i)\}.
$$

For deterministic finite-endpoint samples at $n=25,50,100$, the discrepancy
was 0.08--0.30 EDF units for Uniform and Exponential examples. The asymmetric
Beta$(1,2)$ discrepancy was 0.27--0.33 at $n=50,100$ and 0.84 at $n=25$.
This is adequate for the existing asymptotic EDF convention and does not
justify a gate-derivative correction. Pure-bulk fits retain the old EDF
exactly.

## Fit time

Release builds were timed against `4f3d148`. Bulk-classified fits incur only
sorting and remained within about 3% of baseline for $n=100$--$2000$. Active
experts cost more because the concise implementation directly evaluates the
Gaussian kernels on the 401-point grid:

| active endpoints | $n=100$ | $n=1000$ | $n=2000$ |
|---|---:|---:|---:|
| one | 1.70 ms (4.4x) | 5.93 ms (12.4x) | 10.67 ms (18.5x) |
| two | 2.90 ms (7.5x) | 10.41 ms (23.5x) | 18.84 ms (36.7x) |

The absolute cost remains below 20 ms at $n=2000$. Keep the direct evaluator
for the first PR; optimize only if downstream profiling shows fit time is
material.

## Reproduction

Install the package with a clean native rebuild at each revision and run
`validate_native_boundary_experts.R`, giving each run a distinct method label
and output file. The script writes replication-level ISE/EDF/log-likelihood and
a companion bias/variance summary. Join the two raw files by support,
scenario, sample size, and replication for paired ratios and uncertainty.
