#!/usr/bin/env bash
set -euo pipefail

OUTDIR="${1:-/data/${USER}/gwas_refs/1000g_phase3_allchr}"
BASE_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
MANIFEST="${OUTDIR}/checksums.sha256"
ENV_NAME="${GWAS_PLINK_ENV:-gwas_plink19}"
MAX_JOBS="${MAX_JOBS:-${SLURM_CPUS_PER_TASK:-4}}"
CHR_MEM_MB="${CHR_MEM_MB:-8000}"
MERGE_MEM_MB="${MERGE_MEM_MB:-56000}"
MERGE_THREADS="${MERGE_THREADS:-4}"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"
mkdir -p vcf_raw

have_cmd() {
    command -v "$1" >/dev/null 2>&1
}

if ! have_cmd sha256sum; then
    echo "ERROR: sha256sum not found"
    exit 1
fi

if ! have_cmd wget; then
    echo "ERROR: wget not found"
    exit 1
fi

if ! have_cmd conda; then
    echo "ERROR: conda not found"
    exit 1
fi

init_conda() {
    eval "$(conda shell.bash hook)"
}

ensure_plink_env() {
    init_conda

    if ! conda env list | awk '{print $1}' | grep -Fxq "${ENV_NAME}"; then
        echo "=== Creating conda env: ${ENV_NAME} ==="
        conda create -y -n "${ENV_NAME}" -c bioconda -c conda-forge \
            plink=1.90b6.21 bcftools
    fi

    conda activate "${ENV_NAME}"

    if ! command -v plink >/dev/null 2>&1; then
        echo "=== Installing plink into existing env: ${ENV_NAME} ==="
        conda install -y -n "${ENV_NAME}" -c bioconda -c conda-forge plink=1.90b6.21
    fi

    if ! command -v bcftools >/dev/null 2>&1; then
        echo "=== Installing bcftools into existing env: ${ENV_NAME} ==="
        conda install -y -n "${ENV_NAME}" -c bioconda -c conda-forge bcftools
    fi

    conda activate "${ENV_NAME}"

    command -v plink >/dev/null 2>&1 || {
        echo "ERROR: plink not found after activating ${ENV_NAME}"
        exit 1
    }

    command -v bcftools >/dev/null 2>&1 || {
        echo "ERROR: bcftools not found after activating ${ENV_NAME}"
        exit 1
    }
}

record_checksum() {
    local f="$1"
    local rel
    rel="${f#${OUTDIR}/}"
    rel="${rel#./}"

    mkdir -p "$(dirname "${MANIFEST}")"
    touch "${MANIFEST}"

    grep -Fv "  ${rel}" "${MANIFEST}" > "${MANIFEST}.tmp" || true
    sha256sum "${rel}" >> "${MANIFEST}.tmp"
    mv "${MANIFEST}.tmp" "${MANIFEST}"
}

check_checksum() {
    local f="$1"
    local rel
    rel="${f#${OUTDIR}/}"
    rel="${rel#./}"

    [[ -f "${MANIFEST}" ]] || return 1
    grep -F "  ${rel}" "${MANIFEST}" >/dev/null 2>&1 || return 1
    sha256sum -c <(grep -F "  ${rel}" "${MANIFEST}") >/dev/null 2>&1
}

download_if_needed() {
    local url="$1"
    local out="$2"

    if [[ -s "${out}" ]] && check_checksum "${out}"; then
        echo "[OK] ${out} exists and checksum matches; skipping"
        return 0
    fi

    if [[ -e "${out}" ]]; then
        echo "[WARN] ${out} exists but checksum missing/mismatch; re-downloading"
        rm -f "${out}"
    fi

    echo "[GET] ${out}"
    wget -q --show-progress -O "${out}" "${url}"

    if [[ ! -s "${out}" ]]; then
        echo "ERROR: downloaded file is empty: ${out}"
        exit 1
    fi

    record_checksum "${out}"
}

plink_trio_ok() {
    local prefix="$1"
    [[ -s "${prefix}.bed" && -s "${prefix}.bim" && -s "${prefix}.fam" ]] || return 1
    check_checksum "${prefix}.bed" && check_checksum "${prefix}.bim" && check_checksum "${prefix}.fam"
}

record_plink_trio() {
    local prefix="$1"
    record_checksum "${prefix}.bed"
    record_checksum "${prefix}.bim"
    record_checksum "${prefix}.fam"
}

wait_for_slot() {
    while [ "$(jobs -rp | wc -l)" -ge "${MAX_JOBS}" ]; do
        sleep 2
    done
}

cleanup_chr_intermediates() {
    local chr="$1"
    rm -f "vcf_raw/chr${chr}.dedup.vcf.gz" "vcf_raw/chr${chr}.dedup.vcf.gz.csi"
}

cleanup_after_final_success() {
    echo "=== Cleaning files no longer needed after successful final reference build ==="

    rm -f 1kg_all-merge.*
    rm -f merge_list.txt

    rm -f vcf_raw/chr*.dedup.vcf.gz vcf_raw/chr*.dedup.vcf.gz.csi
    rm -f vcf_raw/chr*.vcf.gz vcf_raw/chr*.vcf.gz.tbi

    rm -f chr*_raw.bed chr*_raw.bim chr*_raw.fam
    rm -f chr*_raw.log chr*_raw.nosex chr*_raw.stdout.log

    echo "=== Cleanup complete ==="
}

convert_chr() {
    local chr="$1"
    local prefix="chr${chr}_raw"
    local vcf="vcf_raw/chr${chr}.vcf.gz"
    local dedup_vcf="vcf_raw/chr${chr}.dedup.vcf.gz"

    init_conda
    conda activate "${ENV_NAME}"

    if [[ ! -s "${vcf}" ]]; then
        echo "ERROR: missing input VCF for chr${chr}: ${vcf}"
        exit 1
    fi

    if plink_trio_ok "${prefix}"; then
        echo "[OK] ${prefix}.{bed,bim,fam} present and checksum matches; skipping conversion"
        cleanup_chr_intermediates "${chr}"
        return 0
    fi

    rm -f "${prefix}.stdout.log"
    echo "[BCFTOOLS] Deduplicating exact duplicate records for chr${chr}" | tee -a "${prefix}.stdout.log"

    rm -f "${dedup_vcf}" "${dedup_vcf}.csi"

    bcftools norm -d all -Oz -o "${dedup_vcf}" "${vcf}" >> "${prefix}.stdout.log" 2>&1
    [[ -s "${dedup_vcf}" ]] || {
        echo "ERROR: bcftools dedup failed for chr${chr}; see ${prefix}.stdout.log"
        exit 1
    }

    bcftools index -f "${dedup_vcf}" >> "${prefix}.stdout.log" 2>&1

    echo "[PLINK] Building ${prefix}" | tee -a "${prefix}.stdout.log"
    rm -f \
        "${prefix}.bed" "${prefix}.bim" "${prefix}.fam" \
        "${prefix}.log" "${prefix}.nosex" \
        "${prefix}"-temporary.*

    plink \
        --vcf "${dedup_vcf}" \
        --double-id \
        --make-bed \
        --biallelic-only \
        --snps-only just-acgt \
        --vcf-half-call m \
        --set-missing-var-ids @:#:\$1:\$2 \
        --memory "${CHR_MEM_MB}" \
        --threads 1 \
        --out "${prefix}" \
        >> "${prefix}.stdout.log" 2>&1

    [[ -s "${prefix}.bed" && -s "${prefix}.bim" && -s "${prefix}.fam" ]] || {
        echo "ERROR: PLINK conversion failed for ${prefix}; see ${prefix}.log and ${prefix}.stdout.log"
        exit 1
    }

    if awk '$2=="."' "${prefix}.bim" | head -n 1 | grep -q .; then
        echo "ERROR: ${prefix}.bim still contains missing variant IDs '.'"
        exit 1
    fi

    record_plink_trio "${prefix}"
    cleanup_chr_intermediates "${chr}"
}

export OUTDIR BASE_URL MANIFEST ENV_NAME MAX_JOBS CHR_MEM_MB MERGE_MEM_MB MERGE_THREADS
export -f have_cmd init_conda ensure_plink_env record_checksum check_checksum
export -f download_if_needed plink_trio_ok record_plink_trio wait_for_slot convert_chr cleanup_chr_intermediates cleanup_after_final_success

echo "=== Output dir: ${OUTDIR} ==="
ensure_plink_env
echo "=== Using conda env: ${ENV_NAME} ==="
echo "=== Parallel jobs: ${MAX_JOBS} ==="
echo "=== Per-chr PLINK memory (MB): ${CHR_MEM_MB} ==="
echo "=== Merge PLINK memory (MB): ${MERGE_MEM_MB} ==="
echo "=== Merge PLINK threads: ${MERGE_THREADS} ==="

echo "=== Downloading 1000 Genomes Phase 3 chr1-22 VCFs if needed ==="
pids=()
for CHR in $(seq 1 22); do
    wait_for_slot
    (
        echo "--- chr${CHR} ---"
        download_if_needed \
            "${BASE_URL}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" \
            "vcf_raw/chr${CHR}.vcf.gz"

        download_if_needed \
            "${BASE_URL}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi" \
            "vcf_raw/chr${CHR}.vcf.gz.tbi"
    ) &
    pids+=($!)
done

for CHR in $(seq 1 22); do
    pid="${pids[$((CHR-1))]}"
    wait "${pid}" || {
        echo "ERROR: chromosome download failed for chr${CHR}"
        exit 1
    }
done

echo "=== Downloading sample panel if needed ==="
download_if_needed \
    "${BASE_URL}/integrated_call_samples_v3.20130502.ALL.panel" \
    "integrated_call_samples_v3.panel"

echo "=== Converting VCFs to PLINK in parallel ==="
pids=()
for CHR in $(seq 1 22); do
    wait_for_slot
    (
        convert_chr "${CHR}"
    ) &
    pids+=($!)
done

for CHR in $(seq 1 22); do
    pid="${pids[$((CHR-1))]}"
    wait "${pid}" || {
        echo "ERROR: chromosome conversion failed for chr${CHR}"
        exit 1
    }
done

echo "=== Validating chr1-22 PLINK trios before merge ==="
for CHR in $(seq 1 22); do
    for ext in bed bim fam; do
        f="chr${CHR}_raw.${ext}"
        [[ -s "${f}" ]] || {
            echo "ERROR: missing required file ${f}"
            exit 1
        }
    done
done

echo "=== Merging chr1-22 into 1kg_all if needed ==="
if plink_trio_ok "1kg_all"; then
    echo "[OK] 1kg_all.{bed,bim,fam} present and checksum matches; skipping merge"
else
    rm -f merge_list.txt
    for CHR in $(seq 2 22); do
        echo "chr${CHR}_raw.bed chr${CHR}_raw.bim chr${CHR}_raw.fam" >> merge_list.txt
    done

    rm -f 1kg_all.bed 1kg_all.bim 1kg_all.fam 1kg_all.log 1kg_all.nosex 1kg_all-merge.* 1kg_all.stdout.log

    plink \
        --bfile chr1_raw \
        --merge-list merge_list.txt \
        --make-bed \
        --memory "${MERGE_MEM_MB}" \
        --threads "${MERGE_THREADS}" \
        --out 1kg_all \
        >> 1kg_all.stdout.log 2>&1

    [[ -s "1kg_all.bed" && -s "1kg_all.bim" && -s "1kg_all.fam" ]] || {
        echo "ERROR: merge failed for 1kg_all"
        exit 1
    }

    if awk '$2=="."' 1kg_all.bim | head -n 1 | grep -q .; then
        echo "ERROR: 1kg_all.bim contains missing variant IDs '.'"
        exit 1
    fi

    record_plink_trio "1kg_all"
    cleanup_after_final_success
fi

echo "=== Done ==="
ls -lh 1kg_all.*
wc -l 1kg_all.fam 1kg_all.bim
echo "Checksum manifest: ${MANIFEST}"