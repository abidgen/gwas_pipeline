#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args      <- commandArgs(trailingOnly = TRUE)
phenotype <- args[1]
covfile   <- args[2]
gwas_pval <- as.numeric(args[3])

pheno <- fread(phenotype)
pcs   <- fread(covfile)

if (!all(c("FID", "IID", "PHE") %in% names(pheno))) {
    stop("Phenotype file must have FID, IID, PHE columns")
}
if (!all(c("FID", "IID") %in% names(pcs))) {
    stop("Covariate/PCA file must have FID, IID columns")
}

merged <- merge(
    pcs,
    pheno,
    by = c("FID", "IID"),
    all.x = TRUE,
    all.y = FALSE
)

merged <- merged[!is.na(PHE)]

# phenotype for PLINK: no header, no quotes
write.table(
    merged[, .(FID, IID, PHE)],
    file = "plink_pheno.txt",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
)

# keep only numeric covariates for PLINK association
keep <- c("FID", "IID")
pc_cols <- grep("^PC[0-9]+$", names(merged), value = TRUE)
keep <- c(keep, pc_cols)

if ("age" %in% names(merged)) keep <- c(keep, "age")
if ("sex" %in% names(merged)) keep <- c(keep, "sex")

covar <- merged[, ..keep]

write.table(
    covar,
    file = "all_covariates.txt",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)

cat("Phenotype and covariates written.\n")
cat("Phenotype rows:", nrow(merged), "\n")
cat("Covariate columns:", paste(names(covar), collapse = ", "), "\n")