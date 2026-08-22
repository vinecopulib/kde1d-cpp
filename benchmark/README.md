# Benchmarks

Run the benchmark suite from any directory with:

```sh
benchmark/run_benchmarks.sh
```

The optional first argument selects the build directory and the second selects
the number of samples per benchmark. For example:

```sh
benchmark/run_benchmarks.sh build-benchmark 50
```

Use the same compiler, build type, machine, and sample count when comparing
commits. Save the output on the base revision and candidate revision, for
example with `benchmark/run_benchmarks.sh | tee benchmark.txt`. The runner
prints the current Git revision along with CMake's compiler information.

The benchmark is deliberately not registered with CTest: timings are for
manual comparisons and never make CI pass or fail.

Evaluation workloads use batches of 10, 100, 1,000, and 10,000 points. They
cover interpolation and integration directly as well as the public PDF and CDF
paths. Quantile timings cover bounded and unbounded continuous fits, discrete
and zero-inflated fits, and the simulation path. Grid construction is measured
both alone and together with interpolation, so an evaluation optimization
cannot hide an excessive setup cost.
