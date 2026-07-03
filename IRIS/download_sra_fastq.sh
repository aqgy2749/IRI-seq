#!/bin/bash
set -euo pipefail

OUT_DIR="/p2/zulab/jtian/data/IRISeq/demo/Data/GSE270383_download_by_aria2c/"

while read acc; do
  [ -z "$acc" ] && continue
  curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${acc}&result=read_run&fields=run_accession,fastq_ftp&format=tsv&download=false" \
  | awk -F'\t' 'NR>1 && $2!="-" {gsub(";", "\nftp://", $2); print "ftp://"$2}'
done < SRR_Acc_List.txt > urls.txt

echo "生成的下载链接预览："
head urls.txt

echo "下载链接数量："
wc -l urls.txt

echo "开始下载 FASTQ..."

# aria2c --all-proxy only accepts an HTTP proxy format: [http://]HOST[:PORT].
# Values such as socks5://127.0.0.1:15780 work for curl, but aria2c rejects them.
proxy_value=""
for candidate in "${ALL_PROXY:-}" "${all_proxy:-}" "${HTTP_PROXY:-}" "${http_proxy:-}" "${HTTPS_PROXY:-}" "${https_proxy:-}"; do
  if [[ "${candidate}" =~ ^(http://)?[^[:space:]/]+(:[0-9]+)?/?$ ]]; then
    proxy_value="${candidate}"
    break
  fi
done

if [[ -n "${proxy_value}" ]]; then
  env -u all_proxy -u ALL_PROXY aria2c --all-proxy="${proxy_value}" -x 12 -s 12 -j 4 -c -d "${OUT_DIR}" -i urls.txt
else
  if [[ -n "${ALL_PROXY:-${all_proxy:-}}" ]]; then
    echo "检测到 aria2c 不支持的 all_proxy/ALL_PROXY，已在本次下载中忽略：${ALL_PROXY:-${all_proxy:-}}" >&2
  fi
  aria2c -x 12 -s 12 -j 4 -c -d "${OUT_DIR}" -i urls.txt
fi
