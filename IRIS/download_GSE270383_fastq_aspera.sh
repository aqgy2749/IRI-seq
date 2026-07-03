#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="/p2/zulab/jtian/data/IRISeq/demo/Data/GSE270383"
LOG_PREFIX="[GSE270383 aspera]"

# Tune this if the connection is unstable. Examples: 100m, 300m, 500m, 1g.
ASPERA_RATE="${ASPERA_RATE:-300m}"

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

cat > GSE270383_Hip_demo_fastq.md5 <<'MD5'
db50e6caa46c1aaf07d0202512f7e9b7  SRR29481264_1.fastq.gz
fedc4fa756a97828a4acc8ca23733f1d  SRR29481264_2.fastq.gz
2e255c3c5704b6cd6fdfe366c19910b1  SRR29481263_1.fastq.gz
b82420ca851918e5036d114cbd0ef098  SRR29481263_2.fastq.gz
MD5

if ! command -v ascp >/dev/null 2>&1; then
    echo "${LOG_PREFIX} [ERROR] ascp not found. Activate/install an Aspera conda env first." >&2
    echo "${LOG_PREFIX} Example: conda create -n aspera -c hcc aspera-cli && conda activate aspera" >&2
    exit 1
fi

find_aspera_key() {
    local candidate

    if [[ -n "${ASPERA_KEY:-}" && -f "${ASPERA_KEY}" ]]; then
        printf '%s\n' "${ASPERA_KEY}"
        return 0
    fi

    for candidate in \
        "${CONDA_PREFIX:-}/etc/asperaweb_id_dsa.openssh" \
        "${CONDA_PREFIX:-}/etc/aspera/asperaweb_id_dsa.openssh" \
        "${HOME}/.aspera/connect/etc/asperaweb_id_dsa.openssh" \
        "/opt/aspera/connect/etc/asperaweb_id_dsa.openssh" \
        "/usr/local/aspera/connect/etc/asperaweb_id_dsa.openssh"
    do
        if [[ -f "${candidate}" ]]; then
            printf '%s\n' "${candidate}"
            return 0
        fi
    done

    return 1
}

ASPERA_KEY_PATH="$(find_aspera_key || true)"
if [[ -z "${ASPERA_KEY_PATH}" ]]; then
    echo "${LOG_PREFIX} [ERROR] Could not find asperaweb_id_dsa.openssh." >&2
    echo "${LOG_PREFIX} Set it manually, for example:" >&2
    echo "${LOG_PREFIX}   ASPERA_KEY=/path/to/asperaweb_id_dsa.openssh bash $0" >&2
    exit 1
fi

download_one() {
    local file="$1"
    local remote_path="$2"
    local expected_md5="$3"

    if [[ -s "${file}" ]] && echo "${expected_md5}  ${file}" | md5sum -c - >/dev/null 2>&1; then
        echo "${LOG_PREFIX} [OK] ${file} already exists and md5 passed"
        return 0
    fi

    echo "${LOG_PREFIX} [START] ${file}"
    echo "${LOG_PREFIX} [KEY] ${ASPERA_KEY_PATH}"
    echo "${LOG_PREFIX} [RATE] ${ASPERA_RATE}"

    ascp -QT -l "${ASPERA_RATE}" -P33001 -k 1 -i "${ASPERA_KEY_PATH}" \
        "era-fasp@fasp.sra.ebi.ac.uk:${remote_path}" .

    echo "${expected_md5}  ${file}" | md5sum -c -
    echo "${LOG_PREFIX} [DONE] ${file}"
}

download_one \
    "SRR29481264_2.fastq.gz" \
    "/vol1/fastq/SRR294/064/SRR29481264/SRR29481264_2.fastq.gz" \
    "fedc4fa756a97828a4acc8ca23733f1d"

download_one \
    "SRR29481263_1.fastq.gz" \
    "/vol1/fastq/SRR294/063/SRR29481263/SRR29481263_1.fastq.gz" \
    "2e255c3c5704b6cd6fdfe366c19910b1"

download_one \
    "SRR29481263_2.fastq.gz" \
    "/vol1/fastq/SRR294/063/SRR29481263/SRR29481263_2.fastq.gz" \
    "b82420ca851918e5036d114cbd0ef098"

md5sum -c GSE270383_Hip_demo_fastq.md5
