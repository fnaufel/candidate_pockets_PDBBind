#!/usr/bin/env bash
#SBATCH --job-name=biosensia-combine-set
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail

# Submit this script from the repository root. SLURM_SUBMIT_DIR is the
# directory in which sbatch was invoked, not the location of this script.
readonly REPO_DIR="${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
cd "${REPO_DIR}"

readonly CONFIG="config/drugclip-combine-set-library.toml"
readonly SOURCE_ROOT="BioSensIA-DC/external/DrugCLIP/data/pdb/combine_set"
readonly DICTIONARY="BioSensIA-DC/external/DrugCLIP/data/dict_pkt.txt"
readonly WORKERS="${SLURM_CPUS_PER_TASK:-1}"

if [[ ! -f "${CONFIG}" || ! -d "${SOURCE_ROOT}" || ! -f "${DICTIONARY}" ]]; then
    printf 'Required input is missing under %s\n' "${REPO_DIR}" >&2
    exit 2
fi

if ! command -v uv >/dev/null 2>&1; then
    printf 'uv is not available in the batch environment PATH\n' >&2
    exit 2
fi

# The environment must be synchronized on a networked node before submission.
# These settings prevent dependency and application-level network access.
export UV_OFFLINE=1
export PYTHONUNBUFFERED=1
export TOKENIZERS_PARALLELISM=false

# Record processing is parallelized by the application. Prevent numerical
# libraries from creating additional thread pools for every worker.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

printf 'Started: %s\n' "$(date --iso-8601=seconds)"
printf 'Job ID: %s\n' "${SLURM_JOB_ID:-unknown}"
printf 'Node: %s\n' "$(hostname)"
printf 'Repository: %s\n' "${REPO_DIR}"
printf 'Git commit: %s\n' "$(git rev-parse HEAD)"
printf 'Workers: %s\n' "${WORKERS}"
uv --version
df -h . data/processed 2>/dev/null || df -h .

time_command=()
if [[ -x /usr/bin/time ]]; then
    time_command=(/usr/bin/time -v)
fi

# Omitting --limit, --pdb-id, and --pdb-ids-file selects every discovered
# data.pkl bundle. --resume makes repeat submission idempotent for a compatible
# completed run and permits rebuilding a compatible interrupted/failed run.
"${time_command[@]}" uv run --offline --no-sync biosensia-pocket-library \
    build-combine-set \
    --config "${CONFIG}" \
    --offline \
    --workers "${WORKERS}" \
    --resume

printf 'Finished: %s\n' "$(date --iso-8601=seconds)"
