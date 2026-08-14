# Run a reproducible, two-sided E1-neck study with separate floating-point
# srvey and rigorous interval-proof products.

# Typical WSL use:
#   MODE=both \
#   SURVEY_BINARY=build/hill4bp_poincare_capd \
#   PROOF_BINARY=build/hill4bp_validate_capd \
#   DELTA_LEVELS="1e-2 1e-3 1e-4 1e-5 1e-6" \
#   SHARD_COUNT=8 JOBS=8 \
#   bash apps/run_professional_neck_study.sh

set -euo pipefail
export LC_ALL=C

launch_dir="${PWD}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(cd "${script_dir}/.." && pwd)"
cd "${project_dir}"

print_help() {
    cat <<'EOF'
Usage: bash apps/run_professional_neck_study.sh [survey|proof|both]

Runs every DELTA_LEVELS value on both the closed and open sides of E1.
The survey is the dense D* CAPD point computation; the proof is the separate
I* CAPD validator. Important environment settings:

  SURVEY_BINARY, PROOF_BINARY     executable paths
  OUTPUT_DIR, STUDY_ID            result location and stable study name
  DELTA_LEVELS                    whitespace-separated positive exact decimals
  NECK_WINDOW_SIGMA, NECK_BRANCH  scaled E1 window (or X_BEGIN and X_END)
  SURVEY_X_COUNT, SURVEY_FRACTION_COUNT, SURVEY_ITERATES
  SURVEY_MAX_ENERGY_ERROR_RATIO    required |H-h|/delta upper limit
  PROOF_X_COUNT, PROOF_FRACTION_COUNT, PROOF_ITERATES
  PROOF_GRID_MODE=points|cover     cover proves an explicitly split rectangle
  PROOF_X_RADIUS, PROOF_FRACTION_RADIUS
  PROOF_REQUIRE_RETURNS=0|1        require every admissible leaf to return
  PROOF_MAX_SUBDIVISION_DEPTH, PROOF_EVENT_TIME_SUBDIVISION_DEPTH
  PROOF_MAX_EVENT_SCAN_NODES, PROOF_MAX_STEP_RETRIES, REQUIRE_COMPLETE=0|1
  SHARD_COUNT, JOBS               deterministic x shards and local concurrency
  SHARD_INDEX                     run one zero-based shard (for a scheduler)
  RESUME=1                        skip matching jobs that have a DONE marker

Each level/side has per-shard logs and metadata. A full-shard invocation also
creates merged CSV tables. Resume settings and executable/configuration hashes
must match the original study.
EOF
}

if [[ "${1:-}" == "help" || "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
    print_help
    exit 0
fi

mode="${1:-${MODE:-both}}"
if [[ $# -gt 1 ]]; then
    echo "Usage: $0 [survey|proof|both]" >&2
    exit 2
fi
case "${mode}" in
    survey|proof|both) ;;
    *) echo "MODE must be survey, proof, or both (got: ${mode})" >&2; exit 2 ;;
esac

resolve_binary() {
    local value="$1"
    case "${value}" in
        /*) printf '%s\n' "${value}" ;;
        *) printf '%s\n' "${launch_dir}/${value}" ;;
    esac
}

survey_binary="$(resolve_binary "${SURVEY_BINARY:-${project_dir}/build/hill4bp_poincare_capd}")"
proof_binary="$(resolve_binary "${PROOF_BINARY:-${project_dir}/build/hill4bp_validate_capd}")"

if [[ "${mode}" != "proof" && ! -x "${survey_binary}" ]]; then
    echo "Survey executable not found or not executable: ${survey_binary}" >&2
    exit 1
fi
if [[ "${mode}" != "survey" && ! -x "${proof_binary}" ]]; then
    echo "Proof executable not found or not executable: ${proof_binary}" >&2
    exit 1
fi

file_sha256() {
    local path="$1"
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "${path}" | awk '{print $1}'
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "${path}" | awk '{print $1}'
    else
        printf '%s\n' unavailable
    fi
}
survey_binary_sha256="not_used"
proof_binary_sha256="not_used"
if [[ "${mode}" != "proof" ]]; then
    survey_binary_sha256="$(file_sha256 "${survey_binary}")"
fi
if [[ "${mode}" != "survey" ]]; then
    proof_binary_sha256="$(file_sha256 "${proof_binary}")"
fi
parameter_file="${project_dir}/Hill_4BP_parameters.cfg"
if [[ ! -f "${parameter_file}" ]]; then
    echo "Shared parameter file not found: ${parameter_file}" >&2
    exit 1
fi
parameter_file_sha256="$(file_sha256 "${parameter_file}")"
runner_sha256="$(file_sha256 "${script_dir}/run_professional_neck_study.sh")"

neck="${NECK:-E1}"
if [[ "${neck}" != "E1" ]]; then
    echo "The rigorous professional workflow currently supports NECK=E1 only." >&2
    exit 2
fi

delta_spec="${DELTA_LEVELS:-1e-2 1e-3 1e-4 1e-5 1e-6}"
read -r -a deltas <<< "${delta_spec}"
if [[ ${#deltas[@]} -eq 0 ]]; then
    echo "DELTA_LEVELS must contain at least one positive magnitude" >&2
    exit 2
fi
for delta in "${deltas[@]}"; do
    if [[ ! "${delta}" =~ ^\+?(0*\.?0*[1-9][0-9]*|[1-9][0-9]*\.?[0-9]*)([eE][+-]?[0-9]+)?$ ]]; then
        echo "Each delta must be a positive, nonzero numeric literal; got: ${delta}" >&2
        exit 2
    fi
done

neck_window_sigma="${NECK_WINDOW_SIGMA:-4}"
neck_branch="${NECK_BRANCH:-positive}"
case "${neck_branch}" in
    positive|negative) ;;
    *) echo "NECK_BRANCH must be positive or negative" >&2; exit 2 ;;
esac

window_args=()
if [[ -n "${X_BEGIN+x}" || -n "${X_END+x}" ]]; then
    if [[ -z "${X_BEGIN+x}" || -z "${X_END+x}" ]]; then
        echo "Set both X_BEGIN and X_END, or neither." >&2
        exit 2
    fi
    window_args=("x_begin=${X_BEGIN}" "x_end=${X_END}")
else
    window_args=(
        "neck_window_sigma=${neck_window_sigma}"
        "neck_branch=${neck_branch}"
    )
fi

survey_x_count="${SURVEY_X_COUNT:-160}"
survey_fraction_count="${SURVEY_FRACTION_COUNT:-61}"
survey_iterates="${SURVEY_ITERATES:-50}"
survey_fraction_lower="${SURVEY_FRACTION_LOWER:--0.95}"
survey_fraction_upper="${SURVEY_FRACTION_UPPER:-0.95}"
survey_solver_order="${SURVEY_SOLVER_ORDER:-24}"
survey_abs_tol="${SURVEY_ABS_TOL:-1e-15}"
survey_rel_tol="${SURVEY_REL_TOL:-1e-13}"
survey_max_step="${SURVEY_MAX_STEP:-0.01}"
survey_svg_stride="${SURVEY_SVG_STRIDE:-10}"
survey_max_energy_error_ratio="${SURVEY_MAX_ENERGY_ERROR_RATIO:-0.01}"

proof_x_count="${PROOF_X_COUNT:-41}"
proof_fraction_count="${PROOF_FRACTION_COUNT:-21}"
proof_iterates="${PROOF_ITERATES:-1}"
proof_fraction_lower="${PROOF_FRACTION_LOWER:--0.9}"
proof_fraction_upper="${PROOF_FRACTION_UPPER:-0.9}"
# Point families are the safe all-purpose default on both sides. A closed
# scaled window contains zero-velocity boundaries, so a continuum cover across
# the entire window necessarily has non-transverse boundary cells. Users opt
# into cover mode on explicitly split, strictly classified rectangles.
proof_grid_mode="${PROOF_GRID_MODE:-points}"
proof_x_radius="${PROOF_X_RADIUS:-0}"
proof_fraction_radius="${PROOF_FRACTION_RADIUS:-0}"
proof_solver_order="${PROOF_SOLVER_ORDER:-30}"
proof_abs_tol="${PROOF_ABS_TOL:-1e-14}"
proof_rel_tol="${PROOF_REL_TOL:-1e-14}"
proof_max_step="${PROOF_MAX_STEP:-0.005}"
proof_max_step_retries="${PROOF_MAX_STEP_RETRIES:-24}"
proof_event_time_subdivision_depth="${PROOF_EVENT_TIME_SUBDIVISION_DEPTH:-30}"
proof_max_event_scan_nodes="${PROOF_MAX_EVENT_SCAN_NODES:-65536}"
proof_max_subdivision_depth="${PROOF_MAX_SUBDIVISION_DEPTH:-12}"
proof_max_leaf_boxes="${PROOF_MAX_LEAF_BOXES:-200000}"
proof_require_returns="${PROOF_REQUIRE_RETURNS:-0}"
proof_progress_every="${PROOF_PROGRESS_EVERY:-25}"
require_complete="${REQUIRE_COMPLETE:-1}"

max_return_time="${MAX_RETURN_TIME:-100}"
max_half_crossings="${MAX_HALF_CROSSINGS:-4}"
y_section="${Y_SECTION:-1e-12}"
collision_radius="${COLLISION_RADIUS:-1e-6}"
outer_radius="${OUTER_RADIUS:-5}"
min_transversality="${MIN_TRANSVERSALITY:-1e-12}"
progress_every="${PROGRESS_EVERY:-10}"
iterate_progress_every="${ITERATE_PROGRESS_EVERY:-10}"

positive_integer() {
    local name="$1" value="$2"
    if [[ ! "${value}" =~ ^[1-9][0-9]*$ ]]; then
        echo "${name} must be a positive integer (got: ${value})" >&2
        exit 2
    fi
}
positive_integer SURVEY_X_COUNT "${survey_x_count}"
positive_integer SURVEY_FRACTION_COUNT "${survey_fraction_count}"
positive_integer SURVEY_ITERATES "${survey_iterates}"
positive_integer PROOF_X_COUNT "${proof_x_count}"
positive_integer PROOF_FRACTION_COUNT "${proof_fraction_count}"
positive_integer PROOF_ITERATES "${proof_iterates}"
positive_integer PROOF_EVENT_TIME_SUBDIVISION_DEPTH "${proof_event_time_subdivision_depth}"
positive_integer PROOF_MAX_EVENT_SCAN_NODES "${proof_max_event_scan_nodes}"
if [[ ! "${proof_max_step_retries}" =~ ^[0-9]+$ ]]; then
    echo "PROOF_MAX_STEP_RETRIES must be a nonnegative integer" >&2
    exit 2
fi
positive_integer PROOF_MAX_LEAF_BOXES "${proof_max_leaf_boxes}"
positive_integer PROOF_PROGRESS_EVERY "${proof_progress_every}"
case "${proof_grid_mode}" in
    cover|points) ;;
    *) echo "PROOF_GRID_MODE must be cover or points" >&2; exit 2 ;;
esac
case "${proof_require_returns}" in
    0|1) ;;
    *) echo "PROOF_REQUIRE_RETURNS must be 0 or 1" >&2; exit 2 ;;
esac
case "${require_complete}" in
    0|1) ;;
    *) echo "REQUIRE_COMPLETE must be 0 or 1" >&2; exit 2 ;;
esac
if [[ ! "${survey_max_energy_error_ratio}" =~ ^\+?(0*\.?0*[1-9][0-9]*|[1-9][0-9]*\.?[0-9]*)([eE][+-]?[0-9]+)?$ ]]; then
    echo "SURVEY_MAX_ENERGY_ERROR_RATIO must be positive and nonzero" >&2
    exit 2
fi

shard_count="${SHARD_COUNT:-1}"
jobs="${JOBS:-1}"
positive_integer SHARD_COUNT "${shard_count}"
positive_integer JOBS "${jobs}"
if [[ "${mode}" != "proof" && ${shard_count} -gt ${survey_x_count} ]]; then
    echo "SHARD_COUNT cannot exceed SURVEY_X_COUNT" >&2
    exit 2
fi
if [[ "${mode}" != "survey" && ${shard_count} -gt ${proof_x_count} ]]; then
    echo "SHARD_COUNT cannot exceed PROOF_X_COUNT" >&2
    exit 2
fi

selected_shard="${SHARD_INDEX:-}"
if [[ -n "${selected_shard}" ]]; then
    if [[ ! "${selected_shard}" =~ ^[0-9]+$ ]] \
        || (( selected_shard >= shard_count )); then
        echo "SHARD_INDEX must be in [0,$((shard_count - 1))]" >&2
        exit 2
    fi
fi

run_id="${STUDY_ID:-$(date -u '+%Y%m%dT%H%M%SZ')_$$}"
output_dir="${OUTPUT_DIR:-${project_dir}/professional_neck_studies/${run_id}}"
resume="${RESUME:-0}"
case "${resume}" in 0|1) ;; *) echo "RESUME must be 0 or 1" >&2; exit 2 ;; esac
manifest="${output_dir}/study.manifest.txt"
if [[ -e "${manifest}" && "${resume}" != "1" ]]; then
    echo "Study already exists: ${output_dir}; set RESUME=1 to continue it." >&2
    exit 1
fi
mkdir -p "${output_dir}"

configuration="${output_dir}/study.config.txt"
requested_configuration="${output_dir}/.study.config.requested.$$"
{
    echo "format=hill4bp_professional_neck_study_config_v1"
    echo "runner_sha256=${runner_sha256}"
    echo "parameter_file=${parameter_file}"
    echo "parameter_file_sha256=${parameter_file_sha256}"
    echo "mode=${mode}"
    echo "neck=${neck}"
    echo "window=${window_args[*]}"
    echo "delta_levels=${delta_spec}"
    echo "shard_count=${shard_count}"
    echo "survey_binary=${survey_binary}"
    echo "survey_binary_sha256=${survey_binary_sha256}"
    echo "survey_x_count=${survey_x_count}"
    echo "survey_fraction_count=${survey_fraction_count}"
    echo "survey_fraction_lower=${survey_fraction_lower}"
    echo "survey_fraction_upper=${survey_fraction_upper}"
    echo "survey_iterates=${survey_iterates}"
    echo "survey_solver_order=${survey_solver_order}"
    echo "survey_abs_tol=${survey_abs_tol}"
    echo "survey_rel_tol=${survey_rel_tol}"
    echo "survey_max_step=${survey_max_step}"
    echo "survey_max_energy_error_ratio=${survey_max_energy_error_ratio}"
    echo "proof_binary=${proof_binary}"
    echo "proof_binary_sha256=${proof_binary_sha256}"
    echo "proof_x_count=${proof_x_count}"
    echo "proof_fraction_count=${proof_fraction_count}"
    echo "proof_fraction_lower=${proof_fraction_lower}"
    echo "proof_fraction_upper=${proof_fraction_upper}"
    echo "proof_grid_mode=${proof_grid_mode}"
    echo "proof_x_radius=${proof_x_radius}"
    echo "proof_fraction_radius=${proof_fraction_radius}"
    echo "proof_iterates=${proof_iterates}"
    echo "proof_solver_order=${proof_solver_order}"
    echo "proof_abs_tol=${proof_abs_tol}"
    echo "proof_rel_tol=${proof_rel_tol}"
    echo "proof_max_step=${proof_max_step}"
    echo "proof_max_step_retries=${proof_max_step_retries}"
    echo "proof_event_time_subdivision_depth=${proof_event_time_subdivision_depth}"
    echo "proof_max_event_scan_nodes=${proof_max_event_scan_nodes}"
    echo "proof_max_subdivision_depth=${proof_max_subdivision_depth}"
    echo "proof_max_leaf_boxes=${proof_max_leaf_boxes}"
    echo "proof_require_returns=${proof_require_returns}"
    echo "require_complete=${require_complete}"
    echo "max_return_time=${max_return_time}"
    echo "max_half_crossings=${max_half_crossings}"
    echo "y_section=${y_section}"
    echo "collision_radius=${collision_radius}"
    echo "outer_radius=${outer_radius}"
    echo "min_transversality=${min_transversality}"
} > "${requested_configuration}"
if [[ -e "${configuration}" ]]; then
    if ! cmp -s "${configuration}" "${requested_configuration}"; then
        rm -f "${requested_configuration}"
        echo "Resume settings do not match ${configuration}; use a new OUTPUT_DIR." >&2
        exit 2
    fi
    rm -f "${requested_configuration}"
else
    mv "${requested_configuration}" "${configuration}"
fi

if [[ ! -e "${manifest}" ]]; then
    manifest_temporary="${output_dir}/.study.manifest.$$"
    {
        echo "format=hill4bp_professional_neck_study_v1"
        echo "runner_sha256=${runner_sha256}"
        echo "parameter_file=${parameter_file}"
        echo "parameter_file_sha256=${parameter_file_sha256}"
        echo "study_id=${run_id}"
        echo "created_utc=$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
        echo "mode=${mode}"
        echo "neck=${neck}"
        echo "neck_branch=${neck_branch}"
        echo "neck_window_sigma=${neck_window_sigma}"
        echo "delta_levels=${delta_spec}"
        echo "shard_count=${shard_count}"
        echo "survey_binary=${survey_binary}"
        echo "survey_binary_sha256=${survey_binary_sha256}"
        echo "proof_binary=${proof_binary}"
        echo "proof_binary_sha256=${proof_binary_sha256}"
        echo "survey_grid=${survey_x_count}x${survey_fraction_count}"
        echo "survey_iterates=${survey_iterates}"
        echo "proof_grid=${proof_x_count}x${proof_fraction_count}"
        echo "proof_grid_mode=${proof_grid_mode}"
        echo "proof_iterates=${proof_iterates}"
        echo "proof_x_radius=${proof_x_radius}"
        echo "proof_fraction_radius=${proof_fraction_radius}"
        echo "require_complete=${require_complete}"
    } > "${manifest_temporary}"
    mv "${manifest_temporary}" "${manifest}"
fi

timeout_program=""
if command -v timeout >/dev/null 2>&1; then
    timeout_program="timeout"
elif command -v gtimeout >/dev/null 2>&1; then
    timeout_program="gtimeout"
fi
preflight_timeout="${PREFLIGHT_TIMEOUT_SECONDS:-120}"
survey_timeout="${SURVEY_TIMEOUT_SECONDS:-86400}"
proof_timeout="${PROOF_TIMEOUT_SECONDS:-86400}"

run_logged() {
    local timeout_seconds="$1" log_path="$2"
    shift 2
    local status
    set +e
    if [[ -n "${timeout_program}" ]]; then
        "${timeout_program}" "${timeout_seconds}" "$@" 2>&1 | tee "${log_path}"
        status=${PIPESTATUS[0]}
    else
        "$@" 2>&1 | tee "${log_path}"
        status=${PIPESTATUS[0]}
    fi
    set -e
    return "${status}"
}

shard_bounds() {
    local count="$1" shard="$2"
    local base=$((count / shard_count))
    local remainder=$((count % shard_count))
    local size start
    if (( shard < remainder )); then
        size=$((base + 1))
        start=$((shard * size))
    else
        size=${base}
        start=$((remainder * (base + 1) + (shard - remainder) * base))
    fi
    printf '%s %s\n' "${start}" "$((start + size - 1))"
}

run_survey_job() {
    local level_label="$1" delta="$2" side="$3" shard="$4"
    local bounds start_index end_index job_dir
    bounds="$(shard_bounds "${survey_x_count}" "${shard}")"
    read -r start_index end_index <<< "${bounds}"
    job_dir="${output_dir}/survey/${level_label}/${side}/shard_$(printf '%04d' "${shard}")"
    if [[ "${resume}" == "1" && -f "${job_dir}/DONE" ]]; then
        if [[ -e "${job_dir}/points.csv" \
            && -s "${job_dir}/returns.csv" \
            && -s "${job_dir}/outcomes.csv" \
            && -s "${job_dir}/metadata.txt" ]]; then
            echo "SKIP completed survey ${level_label} ${side} shard ${shard}"
            return 0
        fi
        echo "REPAIR incomplete survey ${level_label} ${side} shard ${shard}"
    fi
    mkdir -p "${job_dir}"
    rm -f "${job_dir}/DONE" "${job_dir}/FAILED"

    local common_args=(
        "neck=E1" "side=${side}" "delta=${delta}"
        "${window_args[@]}"
        "x_count=${survey_x_count}"
        "start_index=${start_index}" "end_index=${end_index}"
        "velocity_mode=fraction"
        "xdot_lower=${survey_fraction_lower}"
        "xdot_upper=${survey_fraction_upper}"
        "xdot_count=${survey_fraction_count}"
        "iterates=${survey_iterates}" "y_section=${y_section}"
        "solver_order=${survey_solver_order}"
        "abs_tol=${survey_abs_tol}" "rel_tol=${survey_rel_tol}"
        "max_step=${survey_max_step}"
        "max_return_time=${max_return_time}"
        "max_half_crossings=${max_half_crossings}"
        "collision_radius=${collision_radius}" "outer_radius=${outer_radius}"
        "min_transversality=${min_transversality}"
        "progress_every=${progress_every}"
        "iterate_progress_every=${iterate_progress_every}"
        "svg_stride=${survey_svg_stride}"
    )
    echo "PREFLIGHT survey ${level_label} ${side} shard ${shard} [${start_index},${end_index}]"
    if run_logged "${preflight_timeout}" "${job_dir}/preflight.log" \
        "${survey_binary}" "${common_args[@]}" dry_run=1; then
        :
    else
        local status=$?
        printf 'phase=preflight\nexit_status=%s\n' "${status}" \
            > "${job_dir}/FAILED"
        return "${status}"
    fi
    echo "RUN survey ${level_label} ${side} shard ${shard}"
    if run_logged "${survey_timeout}" "${job_dir}/run.log" \
        "${survey_binary}" "${common_args[@]}" \
        "output=${job_dir}/points.csv" \
        "returns=${job_dir}/returns.csv" \
        "outcomes=${job_dir}/outcomes.csv" \
        "zvb=${job_dir}/zero_velocity_boundary.csv" \
        "svg=${job_dir}/preview.svg" \
        "metadata=${job_dir}/metadata.txt"; then
        :
    else
        local status=$?
        printf 'phase=run\nexit_status=%s\n' "${status}" \
            > "${job_dir}/FAILED"
        return "${status}"
    fi
    if [[ ! -e "${job_dir}/points.csv" \
        || ! -s "${job_dir}/returns.csv" \
        || ! -s "${job_dir}/outcomes.csv" \
        || ! -s "${job_dir}/metadata.txt" ]]; then
        echo "missing_required_output" > "${job_dir}/FAILED"
        return 1
    fi
    local energy_ratio
    energy_ratio="$(awk -F= '$1 == "max_abs_energy_error_over_delta" {print $2}' \
        "${job_dir}/metadata.txt" | tail -n 1)"
    if [[ -z "${energy_ratio}" ]]; then
        echo "metadata_missing_energy_error_ratio" > "${job_dir}/FAILED"
        return 1
    fi
    if [[ "${energy_ratio}" != "nan" ]]; then
        if [[ ! "${energy_ratio}" =~ ^\+?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$ ]]; then
            echo "invalid_energy_error_ratio=${energy_ratio}" > "${job_dir}/FAILED"
            return 1
        fi
        if awk -v ratio="${energy_ratio}" -v limit="${survey_max_energy_error_ratio}" \
            'BEGIN { exit !(ratio >= limit) }'; then
            printf 'energy_error_ratio=%s\nrequired_less_than=%s\n' \
                "${energy_ratio}" "${survey_max_energy_error_ratio}" \
                > "${job_dir}/FAILED"
            return 1
        fi
    fi
    printf 'completed_utc=%s\nstart_index=%s\nend_index=%s\n' \
        "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" "${start_index}" "${end_index}" \
        > "${job_dir}/DONE"
}

run_proof_job() {
    local level_label="$1" delta="$2" side="$3" shard="$4"
    local bounds start_index end_index job_dir
    bounds="$(shard_bounds "${proof_x_count}" "${shard}")"
    read -r start_index end_index <<< "${bounds}"
    job_dir="${output_dir}/proof/${level_label}/${side}/shard_$(printf '%04d' "${shard}")"
    if [[ "${resume}" == "1" && -f "${job_dir}/DONE" ]]; then
        if [[ -s "${job_dir}/interval_enclosures.csv" \
            && -s "${job_dir}/proof_metadata.txt" ]]; then
            echo "SKIP completed proof ${level_label} ${side} shard ${shard}"
            return 0
        fi
        echo "REPAIR incomplete proof ${level_label} ${side} shard ${shard}"
    fi
    mkdir -p "${job_dir}"
    rm -f "${job_dir}/DONE" "${job_dir}/FAILED"

    local common_args=(
        "neck=E1" "side=${side}" "delta=${delta}"
        "${window_args[@]}"
        "x_count=${proof_x_count}"
        "start_index=${start_index}" "end_index=${end_index}"
        "grid_mode=${proof_grid_mode}"
        "fraction_lower=${proof_fraction_lower}"
        "fraction_upper=${proof_fraction_upper}"
        "fraction_count=${proof_fraction_count}"
        "x_radius=${proof_x_radius}"
        "fraction_radius=${proof_fraction_radius}"
        "iterates=${proof_iterates}"
        "solver_order=${proof_solver_order}"
        "abs_tol=${proof_abs_tol}" "rel_tol=${proof_rel_tol}"
        "max_step=${proof_max_step}"
        "max_step_retries=${proof_max_step_retries}"
        "max_return_time=${max_return_time}"
        "max_half_crossings=${max_half_crossings}"
        "event_time_subdivision_depth=${proof_event_time_subdivision_depth}"
        "max_event_scan_nodes=${proof_max_event_scan_nodes}"
        "max_subdivision_depth=${proof_max_subdivision_depth}"
        "max_leaf_boxes=${proof_max_leaf_boxes}"
        "collision_radius=${collision_radius}" "outer_radius=${outer_radius}"
        "min_transversality=${min_transversality}"
        "require_complete=${require_complete}"
        "require_returns=${proof_require_returns}"
        "progress_every=${proof_progress_every}"
    )
    echo "PREFLIGHT proof ${level_label} ${side} shard ${shard} [${start_index},${end_index}]"
    if run_logged "${preflight_timeout}" "${job_dir}/preflight.log" \
        "${proof_binary}" "${common_args[@]}" dry_run=1; then
        :
    else
        local status=$?
        printf 'phase=preflight\nexit_status=%s\n' "${status}" \
            > "${job_dir}/FAILED"
        return "${status}"
    fi
    echo "RUN proof ${level_label} ${side} shard ${shard}"
    if run_logged "${proof_timeout}" "${job_dir}/run.log" \
        "${proof_binary}" "${common_args[@]}" \
        "output=${job_dir}/interval_enclosures.csv" \
        "metadata=${job_dir}/proof_metadata.txt"; then
        :
    else
        local status=$?
        printf 'phase=run\nexit_status=%s\n' "${status}" \
            > "${job_dir}/FAILED"
        return "${status}"
    fi
    if [[ ! -s "${job_dir}/interval_enclosures.csv" \
        || ! -s "${job_dir}/proof_metadata.txt" ]]; then
        echo "missing_required_output" > "${job_dir}/FAILED"
        return 1
    fi
    local required_metadata=(
        "rigorous_interval_arithmetic=1"
        "fast_math=0"
        "configured_max_step_enforced_and_checked=1"
        "first_neck_certified=1"
        "side_certified=1"
        "side=${side}"
        "delta_literal=${delta}"
        "grid_mode=${proof_grid_mode}"
        "event_detection=separate_capped_IOdeSolver_dense_curve_scan"
        "set_representation=C0HOTripletonSet"
    )
    if [[ "${require_complete}" == "1" ]]; then
        required_metadata+=("proof_complete=1")
    fi
    local assertion
    for assertion in "${required_metadata[@]}"; do
        if ! grep -F -x -q "${assertion}" "${job_dir}/proof_metadata.txt"; then
            printf 'missing_proof_assertion=%s\n' "${assertion}" \
                > "${job_dir}/FAILED"
            return 1
        fi
    done
    if ! grep -E -x -q \
        'capd_interval_backend=(NATIVE|FILIB|CXSC)' \
        "${job_dir}/proof_metadata.txt"; then
        echo "missing_supported_capd_interval_backend" > "${job_dir}/FAILED"
        return 1
    fi
    printf 'completed_utc=%s\nstart_index=%s\nend_index=%s\nrequire_complete=%s\n' \
        "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" "${start_index}" "${end_index}" \
        "${require_complete}" > "${job_dir}/DONE"
}

declare -a batch_pids=()
overall_failed=0

wait_batch() {
    local pid batch_failed=0
    if (( ${#batch_pids[@]} == 0 )); then
        return 0
    fi
    for pid in "${batch_pids[@]}"; do
        if ! wait "${pid}"; then
            batch_failed=1
        fi
    done
    batch_pids=()
    if (( batch_failed != 0 )); then
        overall_failed=1
    fi
}

launch_job() {
    "$@" &
    batch_pids+=("$!")
    if (( ${#batch_pids[@]} >= jobs )); then
        wait_batch
    fi
}

for level_index in "${!deltas[@]}"; do
    delta="${deltas[${level_index}]}"
    delta_label="${delta#+}"
    delta_label="${delta_label//./p}"
    delta_label="${delta_label//+/_plus_}"
    delta_label="${delta_label//-/_minus_}"
    level_label="level_$(printf '%03d' "${level_index}")_delta_${delta_label}"
    for side in closed open; do
        for ((shard = 0; shard < shard_count; ++shard)); do
            if [[ -n "${selected_shard}" && "${shard}" != "${selected_shard}" ]]; then
                continue
            fi
            if [[ "${mode}" != "proof" ]]; then
                launch_job run_survey_job "${level_label}" "${delta}" "${side}" "${shard}"
                if (( overall_failed != 0 )); then
                    break 3
                fi
            fi
            if [[ "${mode}" != "survey" ]]; then
                launch_job run_proof_job "${level_label}" "${delta}" "${side}" "${shard}"
                if (( overall_failed != 0 )); then
                    break 3
                fi
            fi
        done
    done
done
wait_batch

if (( overall_failed != 0 )); then
    echo "Study incomplete: at least one job failed." >&2
    echo "Use RESUME=1 only when the binary and configuration are unchanged; after rebuilding, use a new OUTPUT_DIR." >&2
    exit 1
fi

merge_headered_csv() {
    local destination="$1"
    shift
    local temporary="${destination}.tmp.$$"
    awk 'FNR == 1 && NR != 1 { next } { print }' "$@" > "${temporary}"
    mv "${temporary}" "${destination}"
}

# When this invocation owns all shards, produce deterministic merged tables in
# shard-index order. Per-shard files and metadata remain the primary audit
# record; these merged files are conveniences for analysis and plotting.
if [[ -z "${selected_shard}" ]]; then
    for level_index in "${!deltas[@]}"; do
        delta="${deltas[${level_index}]}"
        delta_label="${delta#+}"
        delta_label="${delta_label//./p}"
        delta_label="${delta_label//+/_plus_}"
        delta_label="${delta_label//-/_minus_}"
        level_label="level_$(printf '%03d' "${level_index}")_delta_${delta_label}"
        for side in closed open; do
            if [[ "${mode}" != "proof" ]]; then
                survey_side_dir="${output_dir}/survey/${level_label}/${side}"
                points_tmp="${survey_side_dir}/points.csv.tmp.$$"
                : > "${points_tmp}"
                outcome_inputs=()
                return_inputs=()
                for ((shard = 0; shard < shard_count; ++shard)); do
                    shard_dir="${survey_side_dir}/shard_$(printf '%04d' "${shard}")"
                    command cat "${shard_dir}/points.csv" >> "${points_tmp}"
                    outcome_inputs+=("${shard_dir}/outcomes.csv")
                    return_inputs+=("${shard_dir}/returns.csv")
                done
                mv "${points_tmp}" "${survey_side_dir}/points.csv"
                merge_headered_csv "${survey_side_dir}/outcomes.csv" \
                    "${outcome_inputs[@]}"
                merge_headered_csv "${survey_side_dir}/returns.csv" \
                    "${return_inputs[@]}"
                cp "${survey_side_dir}/shard_0000/zero_velocity_boundary.csv" \
                    "${survey_side_dir}/zero_velocity_boundary.csv"
                echo "all_shards_complete=1" > "${survey_side_dir}/COMPLETE"
            fi
            if [[ "${mode}" != "survey" ]]; then
                proof_side_dir="${output_dir}/proof/${level_label}/${side}"
                proof_inputs=()
                for ((shard = 0; shard < shard_count; ++shard)); do
                    shard_dir="${proof_side_dir}/shard_$(printf '%04d' "${shard}")"
                    proof_inputs+=("${shard_dir}/interval_enclosures.csv")
                done
                merge_headered_csv "${proof_side_dir}/interval_enclosures.csv" \
                    "${proof_inputs[@]}"
                {
                    echo "all_shards_complete=1"
                    echo "require_complete=${require_complete}"
                    echo "shard_count=${shard_count}"
                } > "${proof_side_dir}/COMPLETE"
            fi
        done
    done
fi

echo "Study complete: ${output_dir}"
