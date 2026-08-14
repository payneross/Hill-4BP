#!/usr/bin/env bash
# Build both the dense survey and rigorous interval validator against one CAPD
# installation.  Intended for Linux/WSL and usable on macOS as well.

set -euo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(cd "${script_dir}/.." && pwd)"
build_dir="${BUILD_DIR:-${project_dir}/build}"
compiler="${CXX:-c++}"

capd_config="${CAPD_CONFIG:-}"
if [[ -z "${capd_config}" ]]; then
    if command -v capd-config >/dev/null 2>&1; then
        capd_config="$(command -v capd-config)"
    elif [[ -x "${project_dir}/../CAPD/build/bin/capd-config" ]]; then
        capd_config="${project_dir}/../CAPD/build/bin/capd-config"
    else
        echo "CAPD_CONFIG is unset and capd-config was not found." >&2
        echo "Set CAPD_CONFIG=/absolute/path/to/CAPD/build/bin/capd-config" >&2
        exit 1
    fi
fi
if [[ ! -x "${capd_config}" ]]; then
    echo "capd-config is not executable: ${capd_config}" >&2
    exit 1
fi
if ! command -v "${compiler}" >/dev/null 2>&1; then
    echo "C++ compiler not found: ${compiler}" >&2
    exit 1
fi

read -r -a capd_cflags <<< "$("${capd_config}" --cflags)"
read -r -a capd_libs <<< "$("${capd_config}" --libs)"
read -r -a extra_cxxflags <<< "${EXTRA_CXXFLAGS:-}"
common_flags=(
    -std=c++17 -O3 -DNDEBUG
    -frounding-math -fno-fast-math -ffp-contract=off
    -Wall -Wextra -Wpedantic
)

mkdir -p "${build_dir}"
echo "CAPD_CONFIG=${capd_config}"
echo "CXX=${compiler}"
echo "Building dense survey..."
"${compiler}" "${capd_cflags[@]}" "${common_flags[@]}" \
    "${extra_cxxflags[@]}" \
    "${project_dir}/apps/hill4bp_poincare_capd.cpp" \
    "${capd_libs[@]}" -o "${build_dir}/hill4bp_poincare_capd"

echo "Building rigorous interval validator..."
"${compiler}" "${capd_cflags[@]}" "${common_flags[@]}" \
    "${extra_cxxflags[@]}" \
    "${project_dir}/apps/hill4bp_validate_capd.cpp" \
    "${capd_libs[@]}" -o "${build_dir}/hill4bp_validate_capd"

"${build_dir}/hill4bp_poincare_capd" \
    neck=E1 side=open delta=1e-4 neck_window_sigma=1 \
    x_count=1 xdot_count=1 iterates=1 dry_run=1 >/dev/null
"${build_dir}/hill4bp_validate_capd" \
    neck=E1 side=open delta=1e-4 neck_window_sigma=1 \
    grid_mode=points x_count=1 fraction_count=1 iterates=1 \
    x_radius=0 fraction_radius=0 min_transversality=0 \
    dry_run=1 >/dev/null

echo "Built and preflighted:"
echo "  ${build_dir}/hill4bp_poincare_capd"
echo "  ${build_dir}/hill4bp_validate_capd"
