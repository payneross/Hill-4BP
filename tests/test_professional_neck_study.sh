#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
test_dir="$(mktemp -d "${TMPDIR:-/tmp}/hill4bp-professional-study-test.XXXXXX")"
trap 'rm -rf "${test_dir}"' EXIT

bash "${project_dir}/apps/run_professional_neck_study.sh" --help \
    > "${test_dir}/help.txt"
grep -q 'survey|proof|both' "${test_dir}/help.txt"

invocations="${test_dir}/invocations.log"
fake_driver="${test_dir}/fake-driver"
printf '%s\n' \
    '#!/usr/bin/env bash' \
    'set -euo pipefail' \
    'kind="$(basename "$0")"' \
    'mode=run' \
    'output=""' \
    'returns=""' \
    'outcomes=""' \
    'metadata=""' \
    'zvb=""' \
    'svg=""' \
    'for argument in "$@"; do' \
    '    case "${argument}" in' \
    '        dry_run=1) mode=preflight ;;' \
    '        output=*) output="${argument#output=}" ;;' \
    '        returns=*) returns="${argument#returns=}" ;;' \
    '        outcomes=*) outcomes="${argument#outcomes=}" ;;' \
    '        metadata=*) metadata="${argument#metadata=}" ;;' \
    '        zvb=*) zvb="${argument#zvb=}" ;;' \
    '        svg=*) svg="${argument#svg=}" ;;' \
    '    esac' \
    'done' \
    'line="${kind} ${mode}"' \
    'for argument in "$@"; do line="${line} ${argument}"; done' \
    'printf "%s\n" "${line}" >> "${FAKE_INVOCATIONS}"' \
    'if [[ "${mode}" == preflight ]]; then exit 0; fi' \
    'if [[ "${kind}" == fake-survey ]]; then' \
    '    : > "${output}"' \
    '    printf "velocity_row,x_index,return_index\n0,0,0\n" > "${returns}"' \
    '    printf "velocity_row,x_index,status\n0,0,completed\n" > "${outcomes}"' \
    '    printf "model=fake_survey\nmax_abs_energy_error_over_delta=%s\n" "${FAKE_ENERGY_RATIO:-1e-8}" > "${metadata}"' \
    '    : > "${zvb}"' \
    '    : > "${svg}"' \
    'else' \
    '    printf "x_index,fraction_index,iterate,status\n0,0,0,proved\n" > "${output}"' \
    '    side=""; delta=""; grid_mode=""' \
    '    for argument in "$@"; do' \
    '        case "${argument}" in side=*) side="${argument#side=}" ;; delta=*) delta="${argument#delta=}" ;; grid_mode=*) grid_mode="${argument#grid_mode=}" ;; esac' \
    '    done' \
    '    printf "rigorous_interval_arithmetic=1\ncapd_interval_backend=NATIVE\nfast_math=0\nconfigured_max_step_enforced_and_checked=1\nfirst_neck_certified=1\nside_certified=1\nside=%s\ndelta_literal=%s\ngrid_mode=%s\nevent_detection=separate_capped_IOdeSolver_dense_curve_scan\nset_representation=C0HOTripletonSet\nproof_complete=1\n" "${side}" "${delta}" "${grid_mode}" > "${metadata}"' \
    'fi' \
    > "${fake_driver}"
chmod +x "${fake_driver}"
ln -s "${fake_driver}" "${test_dir}/fake-survey"
ln -s "${fake_driver}" "${test_dir}/fake-proof"

study_dir="${test_dir}/study with spaces"
if ! FAKE_INVOCATIONS="${invocations}" \
    OUTPUT_DIR="${study_dir}" STUDY_ID=contract MODE=both \
    SURVEY_BINARY="${test_dir}/fake-survey" \
    PROOF_BINARY="${test_dir}/fake-proof" \
    DELTA_LEVELS="1e-2 1e-4" \
    SURVEY_X_COUNT=5 SURVEY_FRACTION_COUNT=3 SURVEY_ITERATES=2 \
    PROOF_X_COUNT=5 PROOF_FRACTION_COUNT=3 PROOF_ITERATES=1 \
    SHARD_COUNT=2 JOBS=3 \
    bash "${project_dir}/apps/run_professional_neck_study.sh" \
        > "${test_dir}/study.log" 2>&1; then
    cat "${test_dir}/study.log" >&2
    exit 1
fi

[[ "$(grep -c '^fake-survey preflight ' "${invocations}")" -eq 8 ]]
[[ "$(grep -c '^fake-survey run ' "${invocations}")" -eq 8 ]]
[[ "$(grep -c '^fake-proof preflight ' "${invocations}")" -eq 8 ]]
[[ "$(grep -c '^fake-proof run ' "${invocations}")" -eq 8 ]]
[[ "$(grep -c ' side=closed ' "${invocations}")" -eq 16 ]]
[[ "$(grep -c ' side=open ' "${invocations}")" -eq 16 ]]
[[ "$(grep -c ' delta=1e-2 ' "${invocations}")" -eq 16 ]]
[[ "$(grep -c ' delta=1e-4 ' "${invocations}")" -eq 16 ]]
[[ "$(grep -c ' start_index=0 end_index=2 ' "${invocations}")" -eq 16 ]]
[[ "$(grep -c ' start_index=3 end_index=4 ' "${invocations}")" -eq 16 ]]
grep -q ' xdot_count=3 ' "${invocations}"
grep -q ' fraction_count=3 ' "${invocations}"
grep -q ' grid_mode=points ' "${invocations}"
grep -q ' max_event_scan_nodes=65536 ' "${invocations}"
grep -q ' require_complete=1' "${invocations}"
[[ "$(find "${study_dir}/survey" -name DONE -type f | wc -l | tr -d ' ')" -eq 8 ]]
[[ "$(find "${study_dir}/proof" -name DONE -type f | wc -l | tr -d ' ')" -eq 8 ]]
[[ -s "${study_dir}/study.manifest.txt" ]]
[[ "$(find "${study_dir}/survey" -name COMPLETE -type f | wc -l | tr -d ' ')" -eq 4 ]]
[[ "$(find "${study_dir}/proof" -name COMPLETE -type f | wc -l | tr -d ' ')" -eq 4 ]]
merged_survey="${study_dir}/survey/level_000_delta_1e_minus_2/closed/outcomes.csv"
merged_returns="${study_dir}/survey/level_000_delta_1e_minus_2/closed/returns.csv"
merged_proof="${study_dir}/proof/level_000_delta_1e_minus_2/closed/interval_enclosures.csv"
[[ "$(wc -l < "${merged_survey}" | tr -d ' ')" -eq 3 ]]
[[ "$(wc -l < "${merged_returns}" | tr -d ' ')" -eq 3 ]]
[[ "$(wc -l < "${merged_proof}" | tr -d ' ')" -eq 3 ]]

before_resume="$(wc -l < "${invocations}" | tr -d ' ')"
FAKE_INVOCATIONS="${invocations}" \
OUTPUT_DIR="${study_dir}" STUDY_ID=contract MODE=both RESUME=1 \
SURVEY_BINARY="${test_dir}/fake-survey" \
PROOF_BINARY="${test_dir}/fake-proof" \
DELTA_LEVELS="1e-2 1e-4" \
SURVEY_X_COUNT=5 SURVEY_FRACTION_COUNT=3 SURVEY_ITERATES=2 \
PROOF_X_COUNT=5 PROOF_FRACTION_COUNT=3 PROOF_ITERATES=1 \
SHARD_COUNT=2 JOBS=3 \
bash "${project_dir}/apps/run_professional_neck_study.sh" \
    > "${test_dir}/resume.log" 2>&1
after_resume="$(wc -l < "${invocations}" | tr -d ' ')"
[[ "${before_resume}" -eq "${after_resume}" ]]
grep -q '^SKIP completed survey' "${test_dir}/resume.log"
grep -q '^SKIP completed proof' "${test_dir}/resume.log"

if FAKE_INVOCATIONS="${invocations}" \
    OUTPUT_DIR="${study_dir}" STUDY_ID=contract MODE=both RESUME=1 \
    SURVEY_BINARY="${test_dir}/fake-survey" \
    PROOF_BINARY="${test_dir}/fake-proof" \
    DELTA_LEVELS="1e-2 1e-4" \
    SURVEY_X_COUNT=5 SURVEY_FRACTION_COUNT=3 SURVEY_ITERATES=2 \
    PROOF_X_COUNT=5 PROOF_FRACTION_COUNT=3 PROOF_ITERATES=2 \
    SHARD_COUNT=2 JOBS=3 \
    bash "${project_dir}/apps/run_professional_neck_study.sh" \
    > "${test_dir}/mismatch.log" 2>&1; then
    echo "resume with changed settings unexpectedly succeeded" >&2
    exit 1
fi
grep -q 'Resume settings do not match' "${test_dir}/mismatch.log"

repair_target="${study_dir}/proof/level_000_delta_1e_minus_2/closed/shard_0000/proof_metadata.txt"
rm -f "${repair_target}"
before_repair="$(wc -l < "${invocations}" | tr -d ' ')"
FAKE_INVOCATIONS="${invocations}" \
OUTPUT_DIR="${study_dir}" STUDY_ID=contract MODE=both RESUME=1 \
SURVEY_BINARY="${test_dir}/fake-survey" \
PROOF_BINARY="${test_dir}/fake-proof" \
DELTA_LEVELS="1e-2 1e-4" \
SURVEY_X_COUNT=5 SURVEY_FRACTION_COUNT=3 SURVEY_ITERATES=2 \
PROOF_X_COUNT=5 PROOF_FRACTION_COUNT=3 PROOF_ITERATES=1 \
SHARD_COUNT=2 JOBS=3 \
bash "${project_dir}/apps/run_professional_neck_study.sh" \
    > "${test_dir}/repair.log" 2>&1
after_repair="$(wc -l < "${invocations}" | tr -d ' ')"
[[ $((after_repair - before_repair)) -eq 2 ]]
grep -q '^REPAIR incomplete proof' "${test_dir}/repair.log"
[[ -s "${repair_target}" ]]

failing_proof="${test_dir}/failing-proof"
printf '%s\n' \
    '#!/usr/bin/env bash' \
    'for argument in "$@"; do [[ "${argument}" == dry_run=1 ]] && exit 0; done' \
    'exit 7' \
    > "${failing_proof}"
chmod +x "${failing_proof}"
if OUTPUT_DIR="${test_dir}/failed-study" MODE=proof \
    PROOF_BINARY="${failing_proof}" DELTA_LEVELS=1e-3 \
    PROOF_X_COUNT=2 PROOF_FRACTION_COUNT=1 PROOF_ITERATES=1 \
    SHARD_COUNT=2 JOBS=2 \
    bash "${project_dir}/apps/run_professional_neck_study.sh" \
    > "${test_dir}/failed.log" 2>&1; then
    echo "proof failure unexpectedly succeeded" >&2
    exit 1
fi
[[ "$(find "${test_dir}/failed-study/proof" -name FAILED -type f | wc -l | tr -d ' ')" -eq 2 ]]
[[ ! -d "${test_dir}/failed-study/proof/level_000_delta_1e_minus_3/open" ]]
grep -R -q '^exit_status=7$' "${test_dir}/failed-study/proof"
grep -q 'Study incomplete' "${test_dir}/failed.log"

if FAKE_INVOCATIONS="${invocations}" FAKE_ENERGY_RATIO=0.02 \
    OUTPUT_DIR="${test_dir}/unresolved-survey" MODE=survey \
    SURVEY_BINARY="${test_dir}/fake-survey" DELTA_LEVELS=1e-3 \
    SURVEY_X_COUNT=1 SURVEY_FRACTION_COUNT=1 SURVEY_ITERATES=1 \
    SHARD_COUNT=1 JOBS=1 \
    bash "${project_dir}/apps/run_professional_neck_study.sh" \
    > "${test_dir}/unresolved-survey.log" 2>&1; then
    echo "survey with excessive energy error unexpectedly succeeded" >&2
    exit 1
fi
grep -R -q '^energy_error_ratio=0.02$' "${test_dir}/unresolved-survey/survey"

if OUTPUT_DIR="${test_dir}/bad-delta" MODE=survey \
    SURVEY_BINARY="${test_dir}/fake-survey" DELTA_LEVELS=0 \
    bash "${project_dir}/apps/run_professional_neck_study.sh" \
    > "${test_dir}/bad-delta.log" 2>&1; then
    echo "zero delta unexpectedly succeeded" >&2
    exit 1
fi
grep -q 'positive, nonzero' "${test_dir}/bad-delta.log"

echo "Professional neck-study workflow contract checks passed"
