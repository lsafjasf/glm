#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
build_dir="${GLM_PRECISION_BUILD_DIR:-${root_dir}/build-precision-regression}"
configuration="${GLM_PRECISION_CONFIG:-Release}"

cmake_args=(
	-S "${root_dir}"
	-B "${build_dir}"
	-DGLM_BUILD_TESTS=ON
	-DGLM_TEST_ENABLE=ON
	-DCMAKE_BUILD_TYPE="${configuration}"
)

if [[ -n "${CMAKE_GENERATOR:-}" ]]; then
	cmake_args+=(-G "${CMAKE_GENERATOR}")
fi

cmake "${cmake_args[@]}"
cmake --build "${build_dir}" --config "${configuration}" \
	--target test-gtx_precision_regression --parallel
ctest --test-dir "${build_dir}" --build-config "${configuration}" \
	-R '^test-gtx_precision_regression$' --output-on-failure
