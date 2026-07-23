#!/usr/bin/env bash
#SBATCH --job-name=hippo_sra_download
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00
#SBATCH --output=/p1/zulab_users/jtian/my_jupyter/IRI-seq/IRIS/hippo_sra_download_%j.out
#SBATCH --error=/p1/zulab_users/jtian/my_jupyter/IRI-seq/IRIS/hippo_sra_download_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=25211510015@m.fudan.edu.cn

set -Eeuo pipefail

OUT_DIR="/p2/zulab/jtian/data/IRISeq/demo_hippo/data"
SRA_TOOLKIT_URL="https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz"
SRA_TOOLKIT_ROOT="${HOME}/.local/opt/sratoolkit-current"

TOTAL_THREADS="${SLURM_CPUS_PER_TASK:-32}"
N_RUNS=2
THREADS_PER_RUN=$(( TOTAL_THREADS / N_RUNS ))
PIGZ_THREADS_PER_FILE=$(( THREADS_PER_RUN / 2 ))

if (( THREADS_PER_RUN < 2 )); then
    echo "ERROR: at least 4 CPUs are required." >&2
    exit 1
fi

mkdir -p "${OUT_DIR}"

# Prefer node-local Slurm scratch. Fall back to /p2, which has much more free
# space than this server's /tmp filesystem.
if [[ -n "${SLURM_TMPDIR:-}" ]]; then
    SCRATCH_ROOT="${SLURM_TMPDIR}/hippo_sra_${SLURM_JOB_ID:-$$}"
else
    SCRATCH_ROOT="${OUT_DIR}/.hippo_sra_tmp_${SLURM_JOB_ID:-$$}"
fi

mkdir -p "${SCRATCH_ROOT}/sra" "${SCRATCH_ROOT}/fastq" "${SCRATCH_ROOT}/fasterq_tmp"

MAIN_BASHPID="${BASHPID}"
cleanup() {
    # Only the top-level shell may remove shared scratch; background workers
    # must leave it intact until both accessions finish.
    if [[ "${BASHPID}" != "${MAIN_BASHPID}" ]]; then
        return
    fi
    if [[ -n "${SCRATCH_ROOT:-}" && -d "${SCRATCH_ROOT}" ]]; then
        rm -rf -- "${SCRATCH_ROOT}"
    fi
}
trap cleanup EXIT INT TERM

ensure_sra_toolkit() {
    if command -v prefetch >/dev/null 2>&1 && command -v fasterq-dump >/dev/null 2>&1; then
        return
    fi

    # Try common cluster module names before installing a user-local copy.
    if command -v module >/dev/null 2>&1; then
        module load sra-tools 2>/dev/null \
            || module load sratoolkit 2>/dev/null \
            || module load SRA-Toolkit 2>/dev/null \
            || true
    fi

    if command -v prefetch >/dev/null 2>&1 && command -v fasterq-dump >/dev/null 2>&1; then
        return
    fi

    mkdir -p "${SRA_TOOLKIT_ROOT}"
    local archive="${SRA_TOOLKIT_ROOT}/sratoolkit.current-ubuntu64.tar.gz"
    local fasterq_path
    # The current NCBI package exposes bin/fasterq-dump as a symbolic link,
    # so do not restrict the search to regular files (-type f).
    fasterq_path="$(find "${SRA_TOOLKIT_ROOT}" -path '*/bin/fasterq-dump' -print -quit)"

    if [[ -z "${fasterq_path}" ]]; then
        echo "SRA Toolkit was not found; installing the current NCBI binary in ${SRA_TOOLKIT_ROOT}"
        if command -v curl >/dev/null 2>&1; then
            curl -fL --retry 5 --retry-delay 10 \
                -o "${archive}.part" "${SRA_TOOLKIT_URL}"
        elif command -v wget >/dev/null 2>&1; then
            wget --tries=5 --waitretry=10 \
                -O "${archive}.part" "${SRA_TOOLKIT_URL}"
        else
            echo "ERROR: neither curl nor wget is available to install SRA Toolkit." >&2
            exit 1
        fi
        mv -f -- "${archive}.part" "${archive}"
        tar -xzf "${archive}" -C "${SRA_TOOLKIT_ROOT}"
        fasterq_path="$(find "${SRA_TOOLKIT_ROOT}" -path '*/bin/fasterq-dump' -print -quit)"
    fi

    if [[ -z "${fasterq_path}" ]]; then
        echo "ERROR: fasterq-dump was not found after SRA Toolkit installation." >&2
        exit 1
    fi

    export PATH="$(dirname "${fasterq_path}"):${PATH}"
}

ensure_sra_toolkit

for tool in prefetch fasterq-dump pigz; do
    if ! command -v "${tool}" >/dev/null 2>&1; then
        echo "ERROR: required command is unavailable: ${tool}" >&2
        exit 1
    fi
done

echo "Job ID: ${SLURM_JOB_ID:-not_running_under_slurm}"
echo "Host: $(hostname)"
echo "Output directory: ${OUT_DIR}"
echo "Scratch directory: ${SCRATCH_ROOT}"
echo "Total CPUs: ${TOTAL_THREADS}"
echo "Threads per accession: ${THREADS_PER_RUN}"
echo "pigz threads per FASTQ: ${PIGZ_THREADS_PER_FILE}"
echo "SRA Toolkit: $(fasterq-dump --version | head -n 1)"
df -h "${OUT_DIR}" "${SCRATCH_ROOT}"

download_and_convert() {
    local accession="$1"
    local output_prefix="$2"
    local final_r1="${OUT_DIR}/${output_prefix}_R1.fastq.gz"
    local final_r2="${OUT_DIR}/${output_prefix}_R2.fastq.gz"
    local accession_dir="${SCRATCH_ROOT}/sra/${accession}"
    local fastq_dir="${SCRATCH_ROOT}/fastq/${accession}"
    local temp_dir="${SCRATCH_ROOT}/fasterq_tmp/${accession}"
    local r1="${fastq_dir}/${accession}_1.fastq"
    local r2="${fastq_dir}/${accession}_2.fastq"

    if [[ -s "${final_r1}" && -s "${final_r2}" ]]; then
        echo "[${accession}] Final R1/R2 files already exist; skipping."
        return
    fi

    mkdir -p "${fastq_dir}" "${temp_dir}"

    echo "[${accession}] Starting prefetch at $(date --iso-8601=seconds)"
    prefetch "${accession}" --max-size u -O "${accession_dir}"

    echo "[${accession}] Starting fasterq-dump at $(date --iso-8601=seconds)"
    (
        cd "${SCRATCH_ROOT}/sra"
        fasterq-dump "${accession}" \
            --split-files \
            --threads "${THREADS_PER_RUN}" \
            --temp "${temp_dir}" \
            --outdir "${fastq_dir}" \
            --progress
    )

    if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
        echo "ERROR: ${accession} did not produce both paired FASTQ files." >&2
        exit 1
    fi

    echo "[${accession}] Compressing R1/R2 with pigz at $(date --iso-8601=seconds)"
    pigz -1 -p "${PIGZ_THREADS_PER_FILE}" "${r1}" &
    local pigz_r1_pid=$!
    pigz -1 -p "${PIGZ_THREADS_PER_FILE}" "${r2}" &
    local pigz_r2_pid=$!
    wait "${pigz_r1_pid}"
    wait "${pigz_r2_pid}"

    mv -f -- "${r1}.gz" "${final_r1}"
    mv -f -- "${r2}.gz" "${final_r2}"

    echo "[${accession}] Completed at $(date --iso-8601=seconds)"
    ls -lh "${final_r1}" "${final_r2}"
}

# Run both accessions concurrently to use network, CPU, and storage bandwidth.
download_and_convert "SRR29481264" "GSM8341603_cDNA" &
pid_cDNA=$!

download_and_convert "SRR29481263" "GSM8341604_interaction" &
pid_interaction=$!

status=0
if ! wait "${pid_cDNA}"; then
    echo "ERROR: SRR29481264 failed." >&2
    status=1
fi
if ! wait "${pid_interaction}"; then
    echo "ERROR: SRR29481263 failed." >&2
    status=1
fi

if (( status != 0 )); then
    exit "${status}"
fi

echo "All downloads completed successfully."
ls -lh \
    "${OUT_DIR}/GSM8341603_cDNA_R1.fastq.gz" \
    "${OUT_DIR}/GSM8341603_cDNA_R2.fastq.gz" \
    "${OUT_DIR}/GSM8341604_interaction_R1.fastq.gz" \
    "${OUT_DIR}/GSM8341604_interaction_R2.fastq.gz"
