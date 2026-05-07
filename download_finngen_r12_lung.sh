#!/usr/bin/env bash
set -euxo pipefail

BASE="https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release"
MANIFEST_URL="https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv"
OUTDIR="${1:-finngen_r12_lung}"
CHECKSUMS_FILE="checksums.sha256"
MAX_JOBS="${MAX_JOBS:-4}"

ENDPOINTS=(
    "finngen_R12_C3_BRONCHUS_LUNG_EXALLC"
    "finngen_R12_C3_NSCLC_ADENO_EXALLC"
    "finngen_R12_C3_LUNG_NONSMALL_EXALLC"
    "finngen_R12_C3_NSCLC_SQUAM_EXALLC"
)

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

have_cmd() {
    command -v "$1" >/dev/null 2>&1
}

for cmd in sha256sum wget awk zcat sort head column mktemp; do
    if ! have_cmd "$cmd"; then
        echo "ERROR: required command not found: $cmd"
        exit 1
    fi
done

record_checksum() {
    local f="$1"
    local rel="${f#./}"

    touch "${CHECKSUMS_FILE}"
    grep -Fv "  ${rel}" "${CHECKSUMS_FILE}" > "${CHECKSUMS_FILE}.tmp" || true
    sha256sum "${rel}" >> "${CHECKSUMS_FILE}.tmp"
    mv "${CHECKSUMS_FILE}.tmp" "${CHECKSUMS_FILE}"
}

check_checksum() {
    local f="$1"
    local rel="${f#./}"

    [[ -f "${CHECKSUMS_FILE}" ]] || return 1
    grep -F "  ${rel}" "${CHECKSUMS_FILE}" >/dev/null 2>&1 || return 1
    sha256sum -c <(grep -F "  ${rel}" "${CHECKSUMS_FILE}") >/dev/null 2>&1
}

manifest_has_file() {
    local fname="$1"
    awk -F'\t' -v f="$fname" '
        NR == 1 { next }
        $6 ~ ("/" f "$") || $7 ~ ("/" f "$") { found=1; exit }
        END { exit(found ? 0 : 1) }
    ' finngen_R12_manifest.tsv
}

check_manifest_presence() {
    local f="$1"
    local fname
    fname="$(basename "$f")"

    if [[ "$fname" == *.gz.tbi ]]; then
        return 0
    fi

    manifest_has_file "$fname"
}

download_if_needed() {
    local url="$1"
    local out="$2"

    if [[ -s "${out}" ]] && check_checksum "${out}" && check_manifest_presence "${out}"; then
        echo "[OK] ${out} exists, checksum matches, and manifest policy passes; skipping"
        return 0
    fi

    if [[ -e "${out}" ]]; then
        echo "[WARN] ${out} exists but checksum/manifest validation failed; re-downloading"
        rm -f "${out}"
    fi

    echo "[GET] ${out}"
    wget -c -q --show-progress "${url}" -O "${out}"

    if [[ ! -s "${out}" ]]; then
        echo "ERROR: downloaded file is empty: ${out}"
        exit 1
    fi

    record_checksum "${out}"

    if ! check_manifest_presence "${out}"; then
        echo "ERROR: $(basename "${out}") is not listed in finngen_R12_manifest.tsv"
        exit 1
    fi
}

wait_for_slot() {
    while [ "$(jobs -rp | wc -l)" -ge "${MAX_JOBS}" ]; do
        sleep 1
    done
}

download_endpoint_pair() {
    local ep="$1"
    echo "  -> ${ep}"
    download_if_needed "${BASE}/${ep}.gz" "${ep}.gz"
    download_if_needed "${BASE}/${ep}.gz.tbi" "${ep}.gz.tbi"
}

export BASE CHECKSUMS_FILE MAX_JOBS
export -f have_cmd record_checksum check_checksum
export -f manifest_has_file check_manifest_presence
export -f download_if_needed wait_for_slot download_endpoint_pair

echo "=== Downloading FinnGen manifest first ==="
if [[ -s finngen_R12_manifest.tsv ]] && check_checksum finngen_R12_manifest.tsv; then
    echo "[OK] finngen_R12_manifest.tsv exists and checksum matches; skipping"
else
    rm -f finngen_R12_manifest.tsv
    wget -c -q --show-progress "${MANIFEST_URL}" -O finngen_R12_manifest.tsv
    [[ -s finngen_R12_manifest.tsv ]] || { echo "ERROR: manifest download failed"; exit 1; }
    record_checksum finngen_R12_manifest.tsv
fi

echo "=== Downloading FinnGen R12 lung cancer endpoints in parallel ==="
for EP in "${ENDPOINTS[@]}"; do
    wait_for_slot
    download_endpoint_pair "${EP}" &
done
wait

echo "=== Done. Files ==="
ls -lh *.gz | awk '{print $5, $9}'

echo "=== Top hits in C3_NSCLC_ADENO_EXALLC (p < 5e-8) ==="
tmp_hits="$(mktemp)"
zcat finngen_R12_C3_NSCLC_ADENO_EXALLC.gz \
    | awk -F'\t' 'NR==1 || ($7!="NA" && $7+0<5e-8)' \
    | sort -t$'\t' -k7,7g | head -6 > "${tmp_hits}" || true

if command -v column >/dev/null 2>&1; then
    column -t < "${tmp_hits}" || cat "${tmp_hits}"
else
    cat "${tmp_hits}"
fi
rm -f "${tmp_hits}"

echo "=== Checksum manifest ==="
pwd
echo "${CHECKSUMS_FILE}"