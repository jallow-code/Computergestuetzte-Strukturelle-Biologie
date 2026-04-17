#!/bin/bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  ./scripts/process_md_trajectory.sh [run_dir] [base_name]

Defaults:
  run_dir    runs/production/md_results
  base_name  orexin1_md

Outputs:
  <base>_centered.xtc
  <base>_centered_fin.xtc
  <base>_centered_fin_start.gro
  <base>_centered_fin_last.gro
  <base>_centered_prot.xtc
  <base>_centered_prot_start.gro

Notes:
  - Centers on the named index group 'Protein' and writes the full 'System'.
  - Uses the named index group 'Protein_NRE_CHL1_POPC' for the reduced trajectory,
    matching the group number used in after2.sh without depending on numeric IDs.
EOF
}

if [ "${1:-}" = "-h" ] || [ "${1:-}" = "--help" ]; then
    usage
    exit 0
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN_DIR="${1:-${ROOT_DIR}/runs/production/md_results}"
BASE="${2:-orexin1_md}"

TPR="${RUN_DIR}/${BASE}.tpr"
XTC="${RUN_DIR}/${BASE}.xtc"
NDX="${RUN_DIR}/${BASE}.ndx"

CENTERED_XTC="${RUN_DIR}/${BASE}_centered.xtc"
CENTERED_FIN_XTC="${RUN_DIR}/${BASE}_centered_fin.xtc"
CENTERED_FIN_START="${RUN_DIR}/${BASE}_centered_fin_start.gro"
CENTERED_FIN_LAST="${RUN_DIR}/${BASE}_centered_fin_last.gro"
CENTERED_COMPLEX_XTC="${RUN_DIR}/${BASE}_centered_prot.xtc"
CENTERED_COMPLEX_START="${RUN_DIR}/${BASE}_centered_prot_start.gro"

for required in "${TPR}" "${XTC}" "${NDX}"; do
    if [ ! -f "${required}" ]; then
        echo "Missing required input: ${required}" >&2
        exit 1
    fi
done

if ! command -v gmx >/dev/null 2>&1; then
    echo "gmx not found in PATH." >&2
    exit 1
fi

LAST_TIME_PS="$(
    gmx check -f "${XTC}" 2>&1 \
        | tr '\r' '\n' \
        | awk '
            /^Step[[:space:]]+/ {
                frames=$2
                dt=$3
            }
            END {
                if (frames > 0 && dt > 0) {
                    printf "%.3f\n", (frames - 1) * dt
                }
            }
        '
)"

if [ -z "${LAST_TIME_PS}" ]; then
    echo "Could not determine the last frame time from ${XTC}" >&2
    exit 1
fi

rm -f \
    "${CENTERED_XTC}" \
    "${CENTERED_FIN_XTC}" \
    "${CENTERED_FIN_START}" \
    "${CENTERED_FIN_LAST}" \
    "${CENTERED_COMPLEX_XTC}" \
    "${CENTERED_COMPLEX_START}"

echo "Step 1/5: remove PBC jumps while centering on Protein"
printf 'Protein\nSystem\n' | gmx trjconv \
    -s "${TPR}" \
    -f "${XTC}" \
    -o "${CENTERED_XTC}" \
    -pbc nojump \
    -center \
    -n "${NDX}"

echo "Step 2/5: rebuild whole molecules in the centered trajectory"
printf 'System\n' | gmx trjconv \
    -s "${TPR}" \
    -f "${CENTERED_XTC}" \
    -o "${CENTERED_FIN_XTC}" \
    -pbc mol \
    -n "${NDX}"

rm -f "${CENTERED_XTC}"

echo "Step 3/5: write first frame structure"
printf 'System\n' | gmx trjconv \
    -s "${TPR}" \
    -f "${CENTERED_FIN_XTC}" \
    -o "${CENTERED_FIN_START}" \
    -dump 0 \
    -n "${NDX}"

echo "Step 4/5: write last frame structure"
printf 'System\n' | gmx trjconv \
    -s "${TPR}" \
    -f "${CENTERED_FIN_XTC}" \
    -o "${CENTERED_FIN_LAST}" \
    -dump "${LAST_TIME_PS}" \
    -n "${NDX}"

echo "Step 5/5: write reduced trajectory and first frame for Protein_NRE_CHL1_POPC"
printf 'Protein_NRE_CHL1_POPC\n' | gmx trjconv \
    -s "${TPR}" \
    -f "${CENTERED_FIN_XTC}" \
    -o "${CENTERED_COMPLEX_XTC}" \
    -n "${NDX}"

printf 'Protein_NRE_CHL1_POPC\n' | gmx trjconv \
    -s "${TPR}" \
    -f "${CENTERED_FIN_XTC}" \
    -o "${CENTERED_COMPLEX_START}" \
    -dump 0 \
    -n "${NDX}"

echo "Trajectory processing complete."
echo "Processed files in: ${RUN_DIR}"
