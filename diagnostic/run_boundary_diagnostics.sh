#!/bin/sh
set -eu

diagnostic_source_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
diagnostic_build_dir=${1:-"$diagnostic_source_dir/build-diagnostic"}

cmake -S "$diagnostic_source_dir" -B "$diagnostic_build_dir" \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=OFF \
  -DBUILD_DIAGNOSTICS=ON 1>&2
cmake --build "$diagnostic_build_dir" \
  --target kde1d-boundary-diagnostic \
  --config Release 1>&2

diagnostic_executable="$diagnostic_build_dir/bin/kde1d-boundary-diagnostic"
if [ ! -x "$diagnostic_executable" ]; then
  diagnostic_executable="$diagnostic_build_dir/bin/Release/kde1d-boundary-diagnostic.exe"
fi

exec "$diagnostic_executable"
