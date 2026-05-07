#!/usr/bin/env bash
# =============================================================================
# simulate_and_run.sh
# Generates simulated GWAS data and runs the full pipeline on it.
# Run from inside gwas-pipeline/ directory.
# =============================================================================
set -euo pipefail

OUTDIR_SIM="simulated_data"
PLINK_PREFIX="simulated_gwas"
RESULTS="/data/${USER}/gwas_results/simulated"

mkdir -p "$OUTDIR_SIM"
cd "$OUTDIR_SIM"

# ── Step 1: Load R and generate simulated data ────────────────────────────────
echo "=== Step 1: Simulating GWAS data ==="
module load R 2>/dev/null || true   # Biowulf R module

Rscript ../scripts/simulate_gwas_data.R \
    --n_cases    500  \
    --n_controls 500  \
    --n_snps     50000 \
    --n_causal   10   \
    --h2         0.10 \
    --out        ${PLINK_PREFIX} \
    --seed       42

echo ""
echo "=== Step 2: Converting PED/MAP to PLINK binary ==="
module load plink 2>/dev/null || true

plink \
    --file  ${PLINK_PREFIX} \
    --make-bed \
    --out   ${PLINK_PREFIX}

echo ""
echo "=== Simulated data summary ==="
echo "SNPs    : $(wc -l < ${PLINK_PREFIX}.bim)"
echo "Samples : $(wc -l < ${PLINK_PREFIX}.fam)"
echo "Cases   : $(awk '$6==2' ${PLINK_PREFIX}.fam | wc -l)"
echo "Controls: $(awk '$6==1' ${PLINK_PREFIX}.fam | wc -l)"

# ── Step 3: Copy phenotype file to data/ ─────────────────────────────────────
echo ""
echo "=== Step 3: Setting up phenotype file ==="
cp ${PLINK_PREFIX}_phenotypes.csv ../data/phenotypes.csv
echo "Phenotype file copied to data/phenotypes.csv"
head -3 ../data/phenotypes.csv

cd ..

# ── Step 4: Run pipeline on simulated data ────────────────────────────────────
echo ""
echo "=== Step 4: Running pipeline on simulated data ==="
module load nextflow 2>/dev/null || true

nextflow run main.nf \
    -profile biowulflocal \
    -resume \
    --mode prs \
    --input_bed  simulated_data/${PLINK_PREFIX} \
    --phenotype  data/phenotypes.csv \
    --sumstats   /data/${USER}/gwas_refs/finngen_r12_lung/finngen_R12_C3_NSCLC_ADENO_EXALLC.gz \
    --pop        EUR \
    --assoc_test logistic \
    --outdir     ${RESULTS}

echo ""
echo "=== Done! Results in: ${RESULTS} ==="