#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="/p2/zulab/jtian/data/IRISeq/demo/Data/GSE270383"
LOG_PREFIX="[GSE270383 curl]"

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

cat > GSE270383_Hip_demo_fastq.md5 <<'MD5'
db50e6caa46c1aaf07d0202512f7e9b7  SRR29481264_1.fastq.gz
fedc4fa756a97828a4acc8ca23733f1d  SRR29481264_2.fastq.gz
2e255c3c5704b6cd6fdfe366c19910b1  SRR29481263_1.fastq.gz
b82420ca851918e5036d114cbd0ef098  SRR29481263_2.fastq.gz
MD5

download_one() {
    local file="$1"
    local url="$2"
    local expected_md5="$3"
    local attempt

    for attempt in 1 2 3; do
        if [[ -s "${file}" ]] && echo "${expected_md5}  ${file}" | md5sum -c - >/dev/null 2>&1; then
            echo "${LOG_PREFIX} [OK] ${file} already exists and md5 passed"
            return 0
        fi

        if [[ -s "${file}" ]]; then
            echo "${LOG_PREFIX} [RESUME] ${file} attempt ${attempt}"
        else
            echo "${LOG_PREFIX} [START] ${file} attempt ${attempt}"
        fi

        curl -L --fail --retry 30 --retry-delay 20 --connect-timeout 60 \
            --speed-time 120 --speed-limit 1024 \
            -C - \
            -o "${file}" \
            "${url}" || true

        if [[ -s "${file}" ]] && echo "${expected_md5}  ${file}" | md5sum -c -; then
            echo "${LOG_PREFIX} [DONE] ${file}"
            return 0
        fi

        if [[ -e "${file}" ]]; then
            local stamp
            stamp="$(date +%Y%m%d_%H%M%S)"
            echo "${LOG_PREFIX} [BAD_MD5] ${file}; preserving as ${file}.bad_md5_${stamp}" >&2
            mv -f "${file}" "${file}.bad_md5_${stamp}"
        fi
    done

    echo "${LOG_PREFIX} [ERROR] ${file} failed after 3 attempts" >&2
    return 1
}

download_one \
    "SRR29481264_1.fastq.gz" \
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR294/064/SRR29481264/SRR29481264_1.fastq.gz" \
    "db50e6caa46c1aaf07d0202512f7e9b7"

download_one \
    "SRR29481264_2.fastq.gz" \
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR294/064/SRR29481264/SRR29481264_2.fastq.gz" \
    "fedc4fa756a97828a4acc8ca23733f1d"

download_one \
    "SRR29481263_1.fastq.gz" \
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR294/063/SRR29481263/SRR29481263_1.fastq.gz" \
    "2e255c3c5704b6cd6fdfe366c19910b1"

download_one \
    "SRR29481263_2.fastq.gz" \
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR294/063/SRR29481263/SRR29481263_2.fastq.gz" \
    "b82420ca851918e5036d114cbd0ef098"

md5sum -c GSE270383_Hip_demo_fastq.md5
