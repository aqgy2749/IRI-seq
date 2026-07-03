#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="/p2/zulab/jtian/data/IRISeq/demo/Data/GSE270383"
N_JOBS="${N_JOBS:-12}"

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
    local size="$4"
    local part_dir=".parts_${file}"
    local tmp_file="${file}.parallel.tmp"
    local chunk_size=$(( (size + N_JOBS - 1) / N_JOBS ))

    if [[ -s "${file}" ]] && echo "${expected_md5}  ${file}" | md5sum -c - >/dev/null 2>&1; then
        echo "[OK] ${file} already exists and md5 passed"
        return 0
    fi

    mkdir -p "${part_dir}"
    echo "[START] ${file} size=${size} jobs=${N_JOBS}"

    local pids=()
    for ((i = 0; i < N_JOBS; i++)); do
        local start=$(( i * chunk_size ))
        local end=$(( start + chunk_size - 1 ))
        if (( start >= size )); then
            break
        fi
        if (( end >= size )); then
            end=$(( size - 1 ))
        fi
        local part="${part_dir}/part_$(printf '%03d' "${i}")"
        local expected_size=$(( end - start + 1 ))
        if [[ -s "${part}" ]] && [[ "$(stat -c '%s' "${part}")" -eq "${expected_size}" ]]; then
            echo "[SKIP] ${file} part ${i} already complete"
            continue
        fi
        (
            curl -L --fail --retry 20 --retry-delay 10 --connect-timeout 60 \
                -r "${start}-${end}" \
                -o "${part}" \
                "${url}"
            actual_size="$(stat -c '%s' "${part}")"
            if [[ "${actual_size}" -ne "${expected_size}" ]]; then
                echo "[ERROR] ${file} part ${i} size ${actual_size}, expected ${expected_size}" >&2
                exit 1
            fi
        ) &
        pids+=("$!")
    done

    local pid
    for pid in "${pids[@]}"; do
        wait "${pid}"
    done

    : > "${tmp_file}"
    for part in "${part_dir}"/part_*; do
        cat "${part}" >> "${tmp_file}"
    done

    local tmp_size
    tmp_size="$(stat -c '%s' "${tmp_file}")"
    if [[ "${tmp_size}" -ne "${size}" ]]; then
        echo "[ERROR] ${tmp_file} size ${tmp_size}, expected ${size}" >&2
        exit 1
    fi

    echo "${expected_md5}  ${tmp_file}" | md5sum -c -
    mv -f "${tmp_file}" "${file}"
    rm -rf "${part_dir}"
    echo "[DONE] ${file}"
}

download_one \
    "SRR29481264_1.fastq.gz" \
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR294/064/SRR29481264/SRR29481264_1.fastq.gz" \
    "db50e6caa46c1aaf07d0202512f7e9b7" \
    8591779957

download_one \
    "SRR29481264_2.fastq.gz" \
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR294/064/SRR29481264/SRR29481264_2.fastq.gz" \
    "fedc4fa756a97828a4acc8ca23733f1d" \
    10868373512

download_one \
    "SRR29481263_1.fastq.gz" \
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR294/063/SRR29481263/SRR29481263_1.fastq.gz" \
    "2e255c3c5704b6cd6fdfe366c19910b1" \
    2784976545

download_one \
    "SRR29481263_2.fastq.gz" \
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR294/063/SRR29481263/SRR29481263_2.fastq.gz" \
    "b82420ca851918e5036d114cbd0ef098" \
    2563487676

md5sum -c GSE270383_Hip_demo_fastq.md5
