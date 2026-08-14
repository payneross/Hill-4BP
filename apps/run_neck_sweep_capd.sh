# Generate a two-sided, near-neck Poincare sweep with the CAPD driver

# Run from any directory:
#   bash apps/run_neck_sweep_capd.sh /path/to/hill4bp_poincare_capd
# or:
#   CAPD_BINARY=/path/to/hill4bp_poincare_capd bash apps/run_neck_sweep_capd.sh



set -euo pipefail

launch_dir="${PWD}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(cd "${script_dir}/.." && pwd)"
cd "${project_dir}"

binary="${1:-${CAPD_BINARY:-${project_dir}/build/hill4bp_poincare_capd}}"
case "${binary}" in
    /*) ;;
    *) binary="${launch_dir}/${binary}" ;;
esac
run_id="$(date -u '+%Y%m%dT%H%M%SZ')_$$"
output_dir="${OUTPUT_DIR:-${project_dir}/near_neck_capd/${run_id}}"

if [[ ! -x "${binary}" ]]; then
    echo "CAPD executable not found or not executable: ${binary}" >&2
    echo "Build apps/hill4bp_poincare_capd.cpp with CAPD first, or set CAPD_BINARY." >&2
    exit 1
fi

# These defaults make a bounded sweep around the positive E1
# neck. The executable computes an x window proportional to sqrt(delta), so
# every energy receives comparable local resolution.
neck="${NECK:-E1}"
x_count="${X_COUNT:-60}"
neck_window_sigma="${NECK_WINDOW_SIGMA:-4}"
neck_branch="${NECK_BRANCH:-positive}"
velocity_mode="${VELOCITY_MODE:-fraction}"
xdot_lower="${XDOT_LOWER:--0.9}"
xdot_upper="${XDOT_UPPER:-0.9}"
xdot_count="${XDOT_COUNT:-11}"
iterates="${ITERATES:-5}"
max_return_time="${MAX_RETURN_TIME:-100}"
max_half_crossings="${MAX_HALF_CROSSINGS:-4}"
collision_radius="${COLLISION_RADIUS:-1e-6}"
outer_radius="${OUTER_RADIUS:-5}"
solver_order="${SOLVER_ORDER:-20}"
abs_tol="${ABS_TOL:-1e-14}"
rel_tol="${REL_TOL:-1e-12}"
max_step="${MAX_STEP:-0.02}"
min_transversality="${MIN_TRANSVERSALITY:-1e-12}"
y_section="${Y_SECTION:-1e-12}"
progress_every="${PROGRESS_EVERY:-10}"
iterate_progress_every="${ITERATE_PROGRESS_EVERY:-5}"
svg_stride="${SVG_STRIDE:-1}"
run_timeout_seconds="${RUN_TIMEOUT_SECONDS:-1800}"

case "${neck}" in
    E1|E3) ;;
    *)
        echo "NECK must be E1 or E3 (got: ${neck})" >&2
        exit 2
        ;;
esac

grid_args=()
if [[ -n "${X_BEGIN+x}" || -n "${X_END+x}" ]]; then
    if [[ -z "${X_BEGIN+x}" || -z "${X_END+x}" ]]; then
        echo "Set both X_BEGIN and X_END, or neither." >&2
        exit 2
    fi
    grid_args=("x_begin=${X_BEGIN}" "x_end=${X_END}")
elif [[ "${neck}" == "E1" ]]; then
    grid_args=(
        "neck_window_sigma=${neck_window_sigma}"
        "neck_branch=${neck_branch}"
    )
else
    # E3/E4 lie on the y axis, so the E1 x-axis neck scaling is inapplicable.
    echo "NECK=E3 requires explicit X_BEGIN and X_END for the y=0 section." >&2
    exit 2
fi

mkdir -p "${output_dir}"

# Positive magnitudes spaced roughly logarithmically. DELTA_LEVELS may replace
# this whitespace-separated list, but zero is rejected because a neck
# equilibrium has zero speed and is non-transverse to the return section.
delta_spec="${DELTA_LEVELS:-1e-2 1e-3 1e-4 1e-5 1e-6}"
read -r -a deltas <<< "${delta_spec}"
if [[ ${#deltas[@]} -eq 0 ]]; then
    echo "DELTA_LEVELS must contain at least one positive magnitude" >&2
    exit 2
fi

echo "CAPD binary: ${binary}"
echo "Output directory: ${output_dir}"
echo "Neck: ${neck}; deltas: ${deltas[*]}"
echo "Bounded return: tau<=${max_return_time}, rho_collision=${collision_radius}, rho_outer=${outer_radius}"

timeout_program=""
if command -v timeout >/dev/null 2>&1; then
    timeout_program="timeout"
elif command -v gtimeout >/dev/null 2>&1; then
    timeout_program="gtimeout"
else
    echo "Warning: timeout/gtimeout is unavailable; CAPD's internal tau and radius limits still apply." >&2
fi

for delta in "${deltas[@]}"; do
    if [[ ! "${delta}" =~ ^\+?(0*\.?0*[1-9][0-9]*|[1-9][0-9]*\.?[0-9]*)([eE][+-]?[0-9]+)?$ ]]; then
        echo "Each delta must be a positive, nonzero numeric literal; got: ${delta}" >&2
        exit 2
    fi

    # Keep output names portable while preserving the supplied delta.
    delta_label="${delta#+}"
    delta_label="${delta_label//./p}"
    delta_label="${delta_label//+/_plus_}"

    for side in closed open; do
        label="${neck}_${side}_delta_${delta_label}"
        prefix="${output_dir}/${label}"

        echo "=== ${neck} ${side}: |h - H(neck)| = ${delta} ==="
        common_args=(
            "neck=${neck}" "side=${side}" "delta=${delta}"
            "${grid_args[@]}" "x_count=${x_count}"
            "velocity_mode=${velocity_mode}"
            "xdot_lower=${xdot_lower}" "xdot_upper=${xdot_upper}"
            "xdot_count=${xdot_count}"
            "y_section=${y_section}" "iterates=${iterates}"
            "solver_order=${solver_order}" "abs_tol=${abs_tol}"
            "rel_tol=${rel_tol}" "max_step=${max_step}"
            "max_return_time=${max_return_time}"
            "max_half_crossings=${max_half_crossings}"
            "min_transversality=${min_transversality}"
            "collision_radius=${collision_radius}" "outer_radius=${outer_radius}"
            "progress_every=${progress_every}"
            "iterate_progress_every=${iterate_progress_every}"
            "svg_stride=${svg_stride}"
        )

        # Resolve the energy, scaled grid, and job size before constructing CAPD.
        "${binary}" "${common_args[@]}" dry_run=1 \
            2>&1 | tee "${prefix}.preflight.log"

        if [[ -n "${timeout_program}" ]]; then
            "${timeout_program}" "${run_timeout_seconds}" \
                "${binary}" "${common_args[@]}" \
                output="${prefix}.csv" returns="${prefix}_returns.csv" \
                zvb="${prefix}_zvb.csv" \
                outcomes="${prefix}_outcomes.csv" svg="${prefix}.svg" \
                metadata="${prefix}.metadata.txt" 2>&1
        else
            "${binary}" "${common_args[@]}" \
                output="${prefix}.csv" returns="${prefix}_returns.csv" \
                zvb="${prefix}_zvb.csv" \
                outcomes="${prefix}_outcomes.csv" svg="${prefix}.svg" \
                metadata="${prefix}.metadata.txt" 2>&1
        fi | tee "${prefix}.log"
    done
done

echo "Sweep complete: ${output_dir}"
