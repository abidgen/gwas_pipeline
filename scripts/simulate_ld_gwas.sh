#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# LD-aware GWAS simulation using existing chr22 1000G PLINK set
# Backbone: /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_chr22_test
# Output:
#   simulated_ld_data/sim_ld_gwas.{bed,bim,fam}
#   simulated_ld_data/sim_phenotypes.csv
#   simulated_ld_data/causal_snps.txt
# ============================================================

module load R

BACKBONE_PREFIX="${BACKBONE_PREFIX:-/data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_chr22_test}"
OUTDIR="${OUTDIR:-simulated_ld_data}"
OUT_PREFIX="${OUT_PREFIX:-sim_ld_gwas}"

N_SAMPLES="${N_SAMPLES:-500}"
N_CAUSAL="${N_CAUSAL:-10}"
SEED="${SEED:-42}"

mkdir -p "${OUTDIR}"

have_cmd() {
    command -v "$1" >/dev/null 2>&1
}

if ! have_cmd awk; then
    echo "ERROR: awk not found"
    exit 1
fi

if ! have_cmd python3; then
    echo "ERROR: python3 not found"
    exit 1
fi

if ! have_cmd Rscript; then
    echo "ERROR: Rscript not found"
    exit 1
fi

if ! have_cmd plink; then
    echo "ERROR: plink not found on PATH"
    exit 1
fi

for ext in bed bim fam; do
    [[ -s "${BACKBONE_PREFIX}.${ext}" ]] || {
        echo "ERROR: missing backbone file ${BACKBONE_PREFIX}.${ext}"
        exit 1
    }
done

echo "=== Backbone prefix: ${BACKBONE_PREFIX} ==="
echo "=== Output dir: ${OUTDIR} ==="
echo "=== Requested samples: ${N_SAMPLES} ==="
echo "=== Requested causal SNPs: ${N_CAUSAL} ==="
echo "=== Seed: ${SEED} ==="

BACKBONE_N=$(wc -l < "${BACKBONE_PREFIX}.fam")
BACKBONE_M=$(wc -l < "${BACKBONE_PREFIX}.bim")

echo "=== Backbone samples: ${BACKBONE_N} ==="
echo "=== Backbone variants: ${BACKBONE_M} ==="

if [[ "${N_SAMPLES}" -gt "${BACKBONE_N}" ]]; then
    echo "ERROR: requested N_SAMPLES=${N_SAMPLES} exceeds available backbone samples=${BACKBONE_N}"
    exit 1
fi

if [[ "${N_CAUSAL}" -gt "${BACKBONE_M}" ]]; then
    echo "ERROR: requested N_CAUSAL=${N_CAUSAL} exceeds available variants=${BACKBONE_M}"
    exit 1
fi

# 1) Sample a subset of individuals from the backbone
echo "=== Step 1: Sampling target individuals from chr22 backbone ==="
python3 - <<PYEOF
import random

seed = int("${SEED}")
n = int("${N_SAMPLES}")
fam_path = "${BACKBONE_PREFIX}.fam"
out_keep = "${OUTDIR}/sample_keep.txt"

random.seed(seed)

samples = []
with open(fam_path) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            samples.append((parts[0], parts[1]))

chosen = random.sample(samples, n)

with open(out_keep, "w") as out:
    for fid, iid in chosen:
        out.write(f"{fid} {iid}\n")

print(f"Wrote {len(chosen)} sample IDs to {out_keep}")
PYEOF

plink \
  --bfile "${BACKBONE_PREFIX}" \
  --keep "${OUTDIR}/sample_keep.txt" \
  --make-bed \
  --out "${OUTDIR}/${OUT_PREFIX}_orig"


python3 - <<PYEOF
fam_path = "${OUTDIR}/${OUT_PREFIX}_orig.fam"
out_path = "${OUTDIR}/update_ids.txt"

rows = []
with open(fam_path) as f:
    for i, line in enumerate(f, start=1):
        parts = line.strip().split()
        old_fid, old_iid = parts[0], parts[1]
        new_id = f"SIM{i:04d}"
        rows.append((old_fid, old_iid, new_id, new_id))

with open(out_path, "w") as out:
    for old_fid, old_iid, new_fid, new_iid in rows:
        out.write(f"{old_fid} {old_iid} {new_fid} {new_iid}\n")

print(f"Wrote {len(rows)} ID updates to {out_path}")
PYEOF

plink \
  --bfile "${OUTDIR}/${OUT_PREFIX}_orig" \
  --update-ids "${OUTDIR}/update_ids.txt" \
  --make-bed \
  --out "${OUTDIR}/${OUT_PREFIX}"

rm -f \
  "${OUTDIR}/${OUT_PREFIX}_orig.bed" \
  "${OUTDIR}/${OUT_PREFIX}_orig.bim" \
  "${OUTDIR}/${OUT_PREFIX}_orig.fam" \
  "${OUTDIR}/${OUT_PREFIX}_orig.log" \
  "${OUTDIR}/${OUT_PREFIX}_orig.nosex"

# 2) Select causal SNPs from the sampled dataset
echo "=== Step 2: Selecting causal SNPs ==="
python3 - <<PYEOF
import random

seed = int("${SEED}") + 1
n_causal = int("${N_CAUSAL}")
bim_path = "${OUTDIR}/${OUT_PREFIX}.bim"
out_path = "${OUTDIR}/causal_snps.txt"

random.seed(seed)

snps = []
with open(bim_path) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            snps.append(parts[1])

chosen = random.sample(snps, n_causal)

with open(out_path, "w") as out:
    for snp in chosen:
        out.write(f"{snp}\n")

print(f"Wrote {len(chosen)} causal SNP IDs to {out_path}")
PYEOF

# 3) Compute additive scores for the chosen SNPs
echo "=== Step 3: Extracting causal SNP dosages ==="
plink \
  --bfile "${OUTDIR}/${OUT_PREFIX}" \
  --extract "${OUTDIR}/causal_snps.txt" \
  --recode A \
  --out "${OUTDIR}/causal_matrix"

# 4) Simulate phenotype/covariates in R using real genotype structure
echo "=== Step 4: Simulating phenotype and covariates ==="
Rscript - "${OUTDIR}" "${OUT_PREFIX}" "${SEED}" <<'REOF'
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
prefix <- args[2]
seed <- as.integer(args[3])

set.seed(seed)

fam_file <- file.path(outdir, paste0(prefix, ".fam"))
raw_file <- file.path(outdir, "causal_matrix.raw")
csv_file <- file.path(outdir, "sim_phenotypes.csv")

fam <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
colnames(fam) <- c("FID","IID","PAT","MAT","SEX_OLD","PHE_OLD")

raw <- read.table(raw_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

meta_cols <- c("FID","IID","PAT","MAT","SEX","PHENOTYPE")
geno_cols <- setdiff(colnames(raw), meta_cols)

if (length(geno_cols) == 0) {
  stop("No genotype columns found in causal_matrix.raw")
}

G <- as.matrix(raw[, geno_cols, drop = FALSE])
G[is.na(G)] <- 0

# Random causal effect sizes
beta <- rnorm(ncol(G), mean = 0, sd = 0.35)

# Polygenic score from real chr22 genotypes at selected causal SNPs
g_score <- as.vector(scale(G %*% beta))

n <- nrow(fam)
sex <- sample(c(1L, 2L), n, replace = TRUE, prob = c(0.45, 0.55))
age <- round(rnorm(n, mean = 62, sd = 10))
age[age < 40] <- 40
age[age > 85] <- 85
batch <- sample(c("batch1", "batch2", "batch3"), n, replace = TRUE)

# Liability with moderate signal
linpred <- 0.8 * g_score + 0.015 * (age - 62) + ifelse(sex == 2L, 0.15, 0)
noise <- rnorm(n, 0, 1.0)
liability <- linpred + noise

# Balanced case/control split for stable testing
thr <- median(liability)
phe <- ifelse(liability > thr, 2L, 1L)

out <- data.frame(
  FID = fam$FID,
  IID = fam$IID,
  PHE = phe,
  age = age,
  sex = sex,
  batch = batch,
  stringsAsFactors = FALSE
)

write.csv(out, csv_file, row.names = FALSE, quote = FALSE)

cat("Phenotype CSV written:", csv_file, "\n")
cat("Samples:", nrow(out), "\n")
cat("Cases:", sum(out$PHE == 2L), "\n")
cat("Controls:", sum(out$PHE == 1L), "\n")
REOF

echo "=== Step 5: Summary ==="
ls -lh \
  "${OUTDIR}/${OUT_PREFIX}.bed" \
  "${OUTDIR}/${OUT_PREFIX}.bim" \
  "${OUTDIR}/${OUT_PREFIX}.fam" \
  "${OUTDIR}/sim_phenotypes.csv" \
  "${OUTDIR}/causal_snps.txt"

echo "Variants: $(wc -l < "${OUTDIR}/${OUT_PREFIX}.bim")"
echo "Samples : $(wc -l < "${OUTDIR}/${OUT_PREFIX}.fam")"
echo "Cases   : $(awk -F',' 'NR>1 && $3==2 {n++} END{print n+0}' "${OUTDIR}/sim_phenotypes.csv")"
echo "Controls: $(awk -F',' 'NR>1 && $3==1 {n++} END{print n+0}' "${OUTDIR}/sim_phenotypes.csv")"

echo ""
echo "Next step:"
echo ""
echo "nextflow run main.nf \\"
echo "  -profile biowulflocal,conda \\"
echo "  --mode full \\"
echo "  --input_bed ${PWD}/${OUTDIR}/${OUT_PREFIX} \\"
echo "  --phenotype ${PWD}/${OUTDIR}/sim_phenotypes.csv \\"
echo "  --sumstats /data/${USER}/gwas_refs/finngen_r12_lung/finngen_R12_C3_NSCLC_ADENO_EXALLC.gz \\"
echo "  --local_1kg_prefix /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_chr22_test \\"
echo "  --ancestry_ref_prefix /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_chr22_test \\"
echo "  --assoc_test logistic \\"
echo "  --pop EUR \\"
echo "  --outdir results_sim_ld_chr22"