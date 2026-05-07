#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

het <- fread("qc_stats/het.het")
het[, Fstat := (`N(NM)` - `O(HOM)`) / `N(NM)`]
mu  <- mean(het$Fstat)
sdd <- sd(het$Fstat)
fails <- het[Fstat < (mu - 3*sdd) | Fstat > (mu + 3*sdd), .(FID, IID)]
fwrite(fails, "het_fails.txt", sep=" ", col.names=FALSE)
cat("Heterozygosity outliers removed:", nrow(fails), "\n")