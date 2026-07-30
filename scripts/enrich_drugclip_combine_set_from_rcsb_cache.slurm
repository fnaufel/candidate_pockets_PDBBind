#!/usr/bin/env bash
#SBATCH --job-name=biosensia-enrich-rcsb
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail

# Submit this script from the repository root.
readonly REPO_DIR="${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
cd "${REPO_DIR}"

readonly RUN_ID="dc-combine-v1-4ab66717-20356381-0953cd2f-e91bfa03"
readonly PARENT_RUN="data/processed/drugclip_combine_set/${RUN_ID}"
readonly CACHE_DIR="rcsb-cache"
readonly OUTPUT_ROOT="data/processed/drugclip_combine_set"

if [[ ! -f "${PARENT_RUN}/manifest.json" ]]; then
    printf 'Parent run manifest not found: %s/manifest.json\n' "${PARENT_RUN}" >&2
    exit 2
fi

if [[ ! -d "${CACHE_DIR}/objects/sha256" || ! -d "${CACHE_DIR}/request_index/rcsb_mmcif" ]]; then
    printf 'RCSB cache is missing required objects or request index: %s\n' "${CACHE_DIR}" >&2
    exit 2
fi

if ! command -v uv >/dev/null 2>&1; then
    printf 'uv is not available in the batch environment PATH\n' >&2
    exit 2
fi

# The environment must already be synchronized on the management node. Both uv
# and the enrichment workflow are forced offline on the CPU compute node.
export UV_OFFLINE=1
export PYTHONUNBUFFERED=1
export TOKENIZERS_PARALLELISM=false

# Avoid nested numerical-library thread pools. The enrichment loop itself is
# serial; the remaining allocated CPUs are available to Arrow and filesystem I/O.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

printf 'Started: %s\n' "$(date --iso-8601=seconds)"
printf 'Job ID: %s\n' "${SLURM_JOB_ID:-unknown}"
printf 'Node: %s\n' "$(hostname)"
printf 'Repository: %s\n' "${REPO_DIR}"
printf 'Git commit: %s\n' "$(git rev-parse HEAD)"
printf 'Parent run: %s\n' "${PARENT_RUN}"
printf 'Cache: %s\n' "${CACHE_DIR}"
printf 'Derived output root: %s\n' "${OUTPUT_ROOT}"
uv --version
df -h . "${CACHE_DIR}" "${OUTPUT_ROOT}"

time_command=()
if [[ -x /usr/bin/time ]]; then
    time_command=(/usr/bin/time -v)
fi

# Partial enrichment is intentionally not enabled. Every accepted PDB ID must
# resolve to a verified local cache object before the derived run is created.
"${time_command[@]}" uv run --offline --no-sync biosensia-pocket-library \
    enrich-from-cache \
    --run-dir "${PARENT_RUN}" \
    --cache-dir "${CACHE_DIR}" \
    --output-root "${OUTPUT_ROOT}" \
    --resume

printf 'Finished: %s\n' "$(date --iso-8601=seconds)"
