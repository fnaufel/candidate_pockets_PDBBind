#!/usr/bin/env bash
#SBATCH --job-name=biosensia-plan-rcsb
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail

# Submit this script from the repository root.
readonly REPO_DIR="${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
cd "${REPO_DIR}"

readonly RUN_ID="dc-combine-v1-4ab66717-20356381-0953cd2f-e91bfa03"
readonly RUN_DIR="data/processed/drugclip_combine_set/${RUN_ID}"
readonly REQUEST_DIR="rcsb-requests"
readonly REQUEST_PATH="${REQUEST_DIR}/${RUN_ID}-rcsb-request.tsv"

if [[ ! -f "${RUN_DIR}/manifest.json" ]]; then
    printf 'Completed run manifest not found: %s/manifest.json\n' "${RUN_DIR}" >&2
    exit 2
fi

if [[ -e "${REQUEST_PATH}" ]]; then
    printf 'Request already exists and will not be overwritten: %s\n' "${REQUEST_PATH}" >&2
    exit 2
fi

if ! command -v uv >/dev/null 2>&1; then
    printf 'uv is not available in the batch environment PATH\n' >&2
    exit 2
fi

# Dependency resolution and the application are both forced offline. The
# environment must already be synchronized on a networked management node.
export UV_OFFLINE=1
export PYTHONUNBUFFERED=1

printf 'Started: %s\n' "$(date --iso-8601=seconds)"
printf 'Job ID: %s\n' "${SLURM_JOB_ID:-unknown}"
printf 'Node: %s\n' "$(hostname)"
printf 'Repository: %s\n' "${REPO_DIR}"
printf 'Git commit: %s\n' "$(git rev-parse HEAD)"
printf 'Source run: %s\n' "${RUN_DIR}"
printf 'Request output: %s\n' "${REQUEST_PATH}"
uv --version

uv run --offline --no-sync biosensia-pocket-library plan-rcsb \
    --run-dir "${RUN_DIR}" \
    --output "${REQUEST_PATH}"

printf 'Request SHA-256: '
sha256sum "${REQUEST_PATH}"
printf 'Requested PDB IDs: %s\n' "$(( $(wc -l < "${REQUEST_PATH}") - 2 ))"
printf 'Finished: %s\n' "$(date --iso-8601=seconds)"
