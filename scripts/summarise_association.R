#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args      <- commandArgs(trailingOnly=TRUE)
assoc     <- args[1]
gwas_pval <- as.numeric(args[2])

res    <- fread(assoc)
res    <- res[!is.na(P)]
chisq  <- qchisq(res$P, df=1, lower.tail=FALSE)
lambda <- median(chisq) / qchisq(0.5, df=1, lower.tail=FALSE)

cat(sprintf("Genomic inflation factor lambda: %.4f\n", lambda))
cat(sprintf("GW significant hits (p < %g): %d\n", gwas_pval, sum(res$P < gwas_pval)))

sig <- res[P < gwas_pval][order(P)]
fwrite(sig, "gwas_significant_hits.txt", sep="\t")
cat("Association summary complete.\n")