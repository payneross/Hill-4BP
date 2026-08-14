#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
test_dir="$(mktemp -d "${TMPDIR:-/tmp}/hill4bp-sweep-test.XXXXXX")"
trap 'rm -rf "${test_dir}"' EXIT

fake_driver="${test_dir}/fake-capd"
printf '%s\n' \
    '#!/usr/bin/env bash' \
    'set -euo pipefail' \
    'mode=run' \
    'for argument in "$@"; do' \
    '    [[ "${argument}" == "dry_run=1" ]] && mode=preflight' \
    'done' \
    'printf "%s" "${mode}"' \
    'for argument in "$@"; do printf " %s" "${argument}"; done' \
    'printf "\n"' \
    > "${fake_driver}"
chmod +x "${fake_driver}"

log="${test_dir}/sweep.log"
OUTPUT_DIR="${test_dir}/outputs" \
DELTA_LEVELS="1e-2 1e-4" \
X_COUNT=2 XDOT_COUNT=2 ITERATES=1 RUN_TIMEOUT_SECONDS=10 \
bash "${project_dir}/apps/run_neck_sweep_capd.sh" "${fake_driver}" > "${log}" 2>&1

[[ "$(grep -c '^preflight ' "${log}")" -eq 4 ]]
[[ "$(grep -c '^run ' "${log}")" -eq 4 ]]
[[ "$(grep -c ' side=closed ' "${log}")" -eq 4 ]]
[[ "$(grep -c ' side=open ' "${log}")" -eq 4 ]]
[[ "$(grep -c ' delta=1e-2 ' "${log}")" -eq 4 ]]
[[ "$(grep -c ' delta=1e-4 ' "${log}")" -eq 4 ]]
grep -q ' neck_window_sigma=4 ' "${log}"
grep -q ' collision_radius=1e-6 ' "${log}"
grep -q ' outer_radius=5 ' "${log}"
grep -q ' outcomes=' "${log}"
grep -q ' returns=' "${log}"
if grep -q ' DELTA=' "${log}"; then
    echo "sweep passed an uppercase DELTA key" >&2
    exit 1
fi

if OUTPUT_DIR="${test_dir}/bad" DELTA_LEVELS=0 \
    bash "${project_dir}/apps/run_neck_sweep_capd.sh" "${fake_driver}" \
    > "${test_dir}/bad.log" 2>&1; then
    echo "zero delta unexpectedly succeeded" >&2
    exit 1
fi
grep -q 'positive, nonzero' "${test_dir}/bad.log"

echo "Near-neck sweep contract checks passed"
