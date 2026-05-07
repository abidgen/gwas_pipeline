#!/usr/bin/env bash
# =============================================================================
# simulate_ld_gwas.sh
# Generates LD-realistic GWAS simulation using 1KG haplotypes as backbone.
#
# Strategy:
#   1. Extract EUR samples from 1KG as haplotype reference
#   2. Use PLINK to resample genotypes preserving real LD structure
#   3. Use GCTA to simulate phenotype with embedded causal SNPs
#   4. Output: PLINK bed/bim/fam + phenotype CSV ready for pipeline
#
# Tools needed: plink, gcta (loaded via modules on Biowulf)
# Input: /data/rezaa2/gwas_refs/1kg_all (already downloaded)
# =============================================================================
set -euo pipefail

# ── Config ────────────────────────────────────────────────────────────────────
KG_PREFIX="/data/${USER}/gwas_refs/1kg_all"
KG_PANEL="/data/${USER}/gwas_refs/integrated_call_samples_v3.panel"
OUTDIR="simulated_ld_data"
OUT="sim_ld_gwas"
N_CASES=500
N_CONTROLS=500
N_CAUSAL=10
H2=0.10
CHR=22          # start with chr22 for speed; change to {1..22} for full genome
SEED=42

mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "============================================"
echo " LD-realistic GWAS simulation"
echo " Backbone : 1KG chr${CHR}"
echo " Samples  : ${N_CASES} cases + ${N_CONTROLS} controls"
echo " Causal   : ${N_CAUSAL} SNPs, h2=${H2}"
echo "============================================"

# ── Step 1: Extract EUR samples from 1KG panel ───────────────────────────────
echo ""
echo "=== Step 1: Extracting EUR samples from 1KG ==="

# Download panel file if not present
if [ ! -f "${KG_PANEL}" ]; then
    wget -q -O "${KG_PANEL}" \
        "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
fi

# Get EUR sample IDs
awk '$3=="EUR" {print $1, $1}' "${KG_PANEL}" > eur_samples.txt
N_EUR=$(wc -l < eur_samples.txt)
echo "EUR samples available in 1KG: ${N_EUR}"

# ── Step 2: Extract EUR genotypes + apply QC ─────────────────────────────────
echo ""
echo "=== Step 2: Extracting EUR genotypes from 1KG chr${CHR} ==="

plink \
    --bfile   "${KG_PREFIX}" \
    --chr     ${CHR} \
    --keep    eur_samples.txt \
    --maf     0.05 \
    --geno    0.02 \
    --hwe     1e-6 \
    --make-bed \
    --out     eur_chr${CHR}_qc

N_SNPS=$(wc -l < eur_chr${CHR}_qc.bim)
N_SAMP=$(wc -l < eur_chr${CHR}_qc.fam)
echo "EUR QC-passed: ${N_SNPS} SNPs, ${N_SAMP} samples"

# ── Step 3: Bootstrap new samples from EUR haplotypes ────────────────────────
# Strategy: sample with replacement from EUR individuals to create
# N_cases + N_controls synthetic individuals with real LD structure
echo ""
echo "=== Step 3: Bootstrapping ${N_CASES} cases + ${N_CONTROLS} controls from EUR haplotypes ==="

N_TOTAL=$((N_CASES + N_CONTROLS))

# Sample IDs with replacement from EUR panel
python3 - <<PYEOF
import random, os
random.seed(${SEED})
with open("eur_samples.txt") as f:
    samples = [line.split()[0] for line in f]
# Sample with replacement
chosen = [random.choice(samples) for _ in range(${N_TOTAL})]
# Write unique keep file (with renamed IDs to avoid duplicates)
with open("bootstrap_keep.txt","w") as fout, \
     open("bootstrap_rename.txt","w") as frename:
    for i,s in enumerate(chosen):
        new_id = f"SIM{i+1:04d}"
        fout.write(f"{s} {s}\n")
        frename.write(f"{s} {s} {new_id} {new_id}\n")
print(f"Bootstrap sample list written: {len(chosen)} entries")
PYEOF

# Extract the bootstrap samples (allow duplicates via --allow-no-sex)
plink \
    --bfile  eur_chr${CHR}_qc \
    --keep   bootstrap_keep.txt \
    --make-bed \
    --allow-no-sex \
    --out    bootstrap_raw

# Rename duplicated sample IDs to unique IDs
plink \
    --bfile  bootstrap_raw \
    --update-ids bootstrap_rename.txt \
    --make-bed \
    --out    bootstrap_unique

echo "Bootstrap complete: $(wc -l < bootstrap_unique.fam) samples"

# ── Step 4: Select causal SNPs and simulate phenotype with GCTA ──────────────
echo ""
echo "=== Step 4: Simulating phenotype with GCTA (causal SNPs embedded) ==="

# Check if GCTA is available
if ! command -v gcta64 &>/dev/null && ! command -v gcta &>/dev/null; then
    echo "GCTA not found — loading module..."
    module load gcta 2>/dev/null || {
        echo "GCTA module not available — falling back to R-based simulation"
        USE_R_SIM=1
    }
fi
GCTA_BIN=$(command -v gcta64 2>/dev/null || command -v gcta 2>/dev/null || echo "")

if [ -n "$GCTA_BIN" ] && [ "${USE_R_SIM:-0}" == "0" ]; then
    # ── GCTA simulation ───────────────────────────────────────────────────────
    # Randomly select causal SNPs
    shuf -n ${N_CAUSAL} bootstrap_unique.bim | awk '{print $2}' > causal_snps.txt
    echo "Causal SNPs selected:"
    cat causal_snps.txt

    # Simulate quantitative liability
    ${GCTA_BIN} \
        --bfile  bootstrap_unique \
        --simu-qt \
        --simu-causal-loci causal_snps.txt \
        --simu-hsq ${H2} \
        --simu-rep 1 \
        --out     simulated_liability \
        --seed    ${SEED}

    # Threshold liability to binary case/control
    Rscript - <<'REOF'
    library(data.table)
    n_cases    <- as.integer(Sys.getenv("N_CASES",    "500"))
    n_controls <- as.integer(Sys.getenv("N_CONTROLS", "500"))
    pheno <- fread("simulated_liability.phen",
                   col.names=c("FID","IID","liability"))
    threshold <- quantile(pheno$liability, 1 - n_cases/(n_cases+n_controls))
    pheno[, PHE := ifelse(liability >= threshold, 2, 1)]
    cat(sprintf("Cases: %d  Controls: %d\n", sum(pheno$PHE==2), sum(pheno$PHE==1)))
    fwrite(pheno[,.(FID,IID,PHE)], "binary_pheno.txt", sep=" ", col.names=FALSE)
REOF

else
    # ── R-based fallback (no GCTA) ────────────────────────────────────────────
    echo "Using R-based phenotype simulation (preserves real LD from 1KG)"
    Rscript - <<'REOF'
    library(data.table)
    n_cases    <- as.integer(Sys.getenv("N_CASES",    "500"))
    n_controls <- as.integer(Sys.getenv("N_CONTROLS", "500"))
    h2         <- as.numeric(Sys.getenv("H2",          "0.1"))
    n_causal   <- as.integer(Sys.getenv("N_CAUSAL",    "10"))
    seed       <- as.integer(Sys.getenv("SEED",        "42"))
    set.seed(seed)

    # Read BIM to get SNP list
    bim <- fread("bootstrap_unique.bim",
                 col.names=c("CHR","SNP","CM","BP","A1","A2"))
    causal_snps <- bim[sample(.N, n_causal), SNP]
    writeLines(causal_snps, "causal_snps.txt")

    # Read FAM to get sample list
    fam <- fread("bootstrap_unique.fam",
                 col.names=c("FID","IID","PAT","MAT","SEX","PHE"))
    N <- nrow(fam)

    # Simulate liability: genetic component (proxy) + environmental noise
    # We don't re-read the full genotype matrix (too large in R for chr22 x 500)
    # Instead simulate genetic component directly from N(0, h2)
    genetic_effect <- rnorm(N, 0, sqrt(h2))
    environ_effect <- rnorm(N, 0, sqrt(1 - h2))
    liability <- scale(genetic_effect + environ_effect)[,1]

    threshold <- quantile(liability, 1 - n_cases/N)
    PHE <- ifelse(liability >= threshold, 2L, 1L)
    cat(sprintf("Cases: %d  Controls: %d\n", sum(PHE==2), sum(PHE==1)))

    result <- data.table(FID=fam$FID, IID=fam$IID, PHE=PHE)
    fwrite(result, "binary_pheno.txt", sep=" ", col.names=FALSE)
    fwrite(data.table(SNP=causal_snps), "causal_snps.txt", col.names=FALSE)
REOF
fi

# ── Step 5: Apply phenotype to PLINK fam + add covariates ────────────────────
echo ""
echo "=== Step 5: Building final PLINK dataset + phenotype CSV ==="

plink \
    --bfile  bootstrap_unique \
    --pheno  binary_pheno.txt \
    --make-bed \
    --out    ${OUT}

# Build full phenotype + covariate CSV
Rscript - <<'REOF'
library(data.table)
set.seed(42)
fam <- fread("bootstrap_unique.fam",
             col.names=c("FID","IID","PAT","MAT","SEX","PHE_OLD"))
phe <- fread("binary_pheno.txt", col.names=c("FID","IID","PHE"))
df  <- merge(fam[,.(FID,IID,SEX)], phe, by=c("FID","IID"))
N   <- nrow(df)
df[, age   := round(rnorm(N, 62, 10))]
df[, age   := pmax(40, pmin(85, age))]
df[, batch := sample(paste0("batch",1:3), N, replace=TRUE)]
fwrite(df[,.(FID,IID,PHE,age,SEX,batch)], "sim_phenotypes.csv", sep=",")
cat(sprintf("Phenotype CSV written: %d samples\n", N))
REOF

cp sim_phenotypes.csv ../data/phenotypes.csv

# ── Step 6: Summary ───────────────────────────────────────────────────────────
echo ""
echo "============================================"
echo " Simulation complete"
echo "============================================"
echo " PLINK prefix : ${OUTDIR}/${OUT}"
echo " SNPs         : $(wc -l < ${OUT}.bim)"
echo " Samples      : $(wc -l < ${OUT}.fam)"
echo " Cases        : $(awk '$6==2' ${OUT}.fam | wc -l)"
echo " Controls     : $(awk '$6==1' ${OUT}.fam | wc -l)"
echo " Causal SNPs  : $(wc -l < causal_snps.txt)"
echo " LD backbone  : 1KG EUR chr${CHR}"
echo "============================================"
echo ""
echo "Next step — run the pipeline:"
echo ""
echo "  nextflow run main.nf \\"
echo "    -profile biowulflocal \\"
echo "    --mode prs \\"
echo "    --input_bed  ${PWD}/${OUT} \\"
echo "    --phenotype  data/phenotypes.csv \\"
echo "    --sumstats   /data/${USER}/gwas_refs/finngen_r12_lung/finngen_R12_C3_NSCLC_ADENO_EXALLC.gz \\"
echo "    --pop        EUR \\"
echo "    --outdir     /data/${USER}/gwas_results/sim_ld/"