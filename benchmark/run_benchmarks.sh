#!/bin/sh
set -eu

benchmark_source_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
benchmark_build_dir=${1:-"$benchmark_source_dir/build-benchmark"}
benchmark_samples=${2:-20}

benchmark_revision=$(git -C "$benchmark_source_dir" rev-parse --short HEAD 2>/dev/null || true)
if [ -n "$benchmark_revision" ]; then
  printf 'kde1d revision: %s\n' "$benchmark_revision"
fi

cmake -S "$benchmark_source_dir" -B "$benchmark_build_dir" \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=OFF \
  -DBUILD_BENCHMARKS=ON
cmake --build "$benchmark_build_dir" --target kde1d-benchmark --config Release

benchmark_executable="$benchmark_build_dir/bin/kde1d-benchmark"
if [ ! -x "$benchmark_executable" ]; then
  benchmark_executable="$benchmark_build_dir/bin/Release/kde1d-benchmark.exe"
fi

exec "$benchmark_executable" '[!benchmark]' \
  --benchmark-samples "$benchmark_samples" \
  --benchmark-resamples 1000
