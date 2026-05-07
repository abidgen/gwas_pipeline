#!/usr/bin/env Rscript
# =============================================================================
# simulate_gwas_data.R
# Simulates realistic GWAS data for pipeline testing:
#   - N samples (cases + controls)
#   - M SNPs across chromosomes 1-22
#   - A few causal SNPs embedded with real effect sizes
#   - Matching phenotype + covariate file
#   - Outputs PLINK PED/MAP → converted to BED/BIM/FAM
# =============================================================================
suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
})

opt_list <- list(
    make_option("--n_cases",    type="integer", default=500,   help="Number of cases"),
    make_option("--n_controls", type="integer", default=500,   help="Number of controls"),
    make_option("--n_snps",     type="integer", default=50000, help="Number of SNPs"),
    make_option("--n_causal",   type="integer", default=10,    help="Number of causal SNPs"),
    make_option("--h2",         type="double",  default=0.1,   help="SNP heritability"),
    make_option("--out",        type="character", default="simulated_gwas", help="Output prefix"),
    make_option("--seed",       type="integer", default=42,    help="Random seed")
)
opt <- parse_args(OptionParser(option_list=opt_list))
set.seed(opt$seed)

N       <- opt$n_cases + opt$n_controls
M       <- opt$n_snps
n_causal <- opt$n_causal
h2      <- opt$h2

cat(sprintf("Simulating %d cases + %d controls, %d SNPs, %d causal, h2=%.2f\n",
    opt$n_cases, opt$n_controls, M, n_causal, h2))

# ── 1. SNP metadata ───────────────────────────────────────────────────────────
# Distribute SNPs across chr1-22 proportional to chromosome length
chr_lengths <- c(249,243,198,191,181,171,159,146,141,136,135,134,
                 115,107,102,90,83,80,59,63,47,51) # approx Mb
chr_weights  <- chr_lengths / sum(chr_lengths)
snp_chr      <- sample(1:22, M, replace=TRUE, prob=chr_weights)
snp_pos      <- sapply(snp_chr, function(c) sample(1:(chr_lengths[c]*1e6), 1))
snp_ids      <- paste0("rs", sample(1e7:9e7, M, replace=FALSE))
maf          <- runif(M, 0.05, 0.45)   # minor allele frequencies
alleles      <- matrix(sample(c("A","T","C","G"), M*2, replace=TRUE), nrow=M)
# ensure ref != alt
same <- alleles[,1] == alleles[,2]
alleles[same, 2] <- ifelse(alleles[same,1]=="A", "T",
                    ifelse(alleles[same,1]=="T", "A",
                    ifelse(alleles[same,1]=="C", "G", "C")))

# ── 2. Genotype matrix (N x M), 0/1/2 additive coding ────────────────────────
cat("Generating genotype matrix...\n")
# Each SNP: HWE genotype frequencies
G <- matrix(0L, nrow=N, ncol=M)
for (j in 1:M) {
    p <- maf[j]
    probs <- c((1-p)^2, 2*p*(1-p), p^2)
    G[,j] <- sample(0:2, N, replace=TRUE, prob=probs)
}

# ── 3. Simulate phenotype with causal SNPs ────────────────────────────────────
cat("Simulating phenotype...\n")
causal_idx    <- sample(1:M, n_causal)
causal_betas  <- rnorm(n_causal, mean=0, sd=sqrt(h2/n_causal))

# Genetic liability
liability <- G[, causal_idx] %*% causal_betas
# Add environmental noise
env_var   <- 1 - h2
liability <- liability + rnorm(N, 0, sqrt(env_var))
liability <- scale(liability)[,1]

# Assign cases/controls using liability threshold
# Top n_cases get PHE=2 (cases), rest get PHE=1 (controls)
threshold   <- quantile(liability, 1 - opt$n_cases/N)
PHE         <- ifelse(liability >= threshold, 2, 1)

cat(sprintf("Cases: %d  Controls: %d\n", sum(PHE==2), sum(PHE==1)))

# ── 4. Sample metadata ────────────────────────────────────────────────────────
FID <- paste0("FAM", sprintf("%04d", 1:N))
IID <- paste0("IND", sprintf("%04d", 1:N))
age <- round(rnorm(N, mean=62, sd=10))
age <- pmax(40, pmin(85, age))
sex <- sample(1:2, N, replace=TRUE, prob=c(0.45, 0.55))  # slight female excess (never-smoker lung ca.)
batch <- sample(paste0("batch", 1:3), N, replace=TRUE)

# ── 5. Write PLINK PED file ───────────────────────────────────────────────────
cat("Writing PLINK PED file...\n")
# PED format: FID IID PAT MAT SEX PHE [SNP1_A1 SNP1_A2 SNP2_A1 SNP2_A2 ...]
ped <- data.table(
    FID = FID, IID = IID,
    PAT = 0, MAT = 0,
    SEX = sex,
    PHE = PHE
)

# Encode genotypes as allele pairs
allele_cols <- list()
for (j in 1:M) {
    a1 <- alleles[j,1]  # ref
    a2 <- alleles[j,2]  # alt
    geno_str <- ifelse(G[,j]==0, paste(a1,a1),
                ifelse(G[,j]==1, paste(a1,a2),
                                 paste(a2,a2)))
    allele_cols[[j]] <- geno_str
}
geno_mat <- do.call(cbind, allele_cols)
ped_full <- cbind(as.data.frame(ped), as.data.frame(geno_mat))
fwrite(ped_full, paste0(opt$out, ".ped"), sep=" ", col.names=FALSE, quote=FALSE)

# ── 6. Write PLINK MAP file ───────────────────────────────────────────────────
cat("Writing PLINK MAP file...\n")
map <- data.table(
    CHR = snp_chr,
    SNP = snp_ids,
    CM  = 0,
    BP  = snp_pos
)
fwrite(map, paste0(opt$out, ".map"), sep="\t", col.names=FALSE, quote=FALSE)

# ── 7. Write phenotype + covariate CSV ───────────────────────────────────────
cat("Writing phenotype file...\n")
pheno <- data.table(FID=FID, IID=IID, PHE=PHE, age=age, sex=sex, batch=batch)
fwrite(pheno, paste0(opt$out, "_phenotypes.csv"), sep=",")

# ── 8. Write causal SNP truth file (for validation) ──────────────────────────
truth <- data.table(
    SNP   = snp_ids[causal_idx],
    CHR   = snp_chr[causal_idx],
    BP    = snp_pos[causal_idx],
    BETA  = causal_betas,
    MAF   = maf[causal_idx]
)
fwrite(truth, paste0(opt$out, "_causal_snps.txt"), sep="\t")

cat("\n=== Simulation complete ===\n")
cat(sprintf("PED file  : %s.ped  (%.1f MB)\n", opt$out, file.size(paste0(opt$out,".ped"))/1e6))
cat(sprintf("MAP file  : %s.map\n", opt$out))
cat(sprintf("Phenotype : %s_phenotypes.csv\n", opt$out))
cat(sprintf("Truth     : %s_causal_snps.txt\n", opt$out))
cat(sprintf("\nConvert to PLINK binary with:\n"))
cat(sprintf("  plink --file %s --make-bed --out %s\n", opt$out, opt$out))
