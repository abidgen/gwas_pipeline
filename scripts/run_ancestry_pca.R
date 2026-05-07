#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(optparse)
})

opt_list <- list(
    make_option("--n_pcs", type="integer", default=10),
    make_option("--pop",   type="character", default="EUR")
)
opt <- parse_args(OptionParser(option_list=opt_list))

evec <- fread("pca/merged_pca.eigenvec",
              col.names=c("FID","IID", paste0("PC", 1:opt$n_pcs)))
panel <- fread("1kg.panel", header=TRUE)
setnames(panel, c("sample","pop","super_pop","gender"))

evec[, ancestry := ifelse(IID %in% panel$sample,
                          panel[match(IID, sample), super_pop],
                          "STUDY")]

# PCA plot
p <- ggplot(evec, aes(PC1, PC2, color=ancestry)) +
    geom_point(alpha=0.6, size=1) +
    theme_bw() +
    labs(title="PCA: study vs 1KG reference", color="Population")
ggsave("pca/pca_plot.png", p, width=8, height=6, dpi=150)

# Ancestry filter: keep study samples within 6 SD of target pop centroid
ref_pop <- evec[ancestry == opt$pop]
mu1 <- mean(ref_pop$PC1); sd1 <- sd(ref_pop$PC1)
mu2 <- mean(ref_pop$PC2); sd2 <- sd(ref_pop$PC2)

study <- evec[ancestry == "STUDY"]
keep  <- study[
    PC1 > mu1 - 6*sd1 & PC1 < mu1 + 6*sd1 &
    PC2 > mu2 - 6*sd2 & PC2 < mu2 + 6*sd2
]
cat("Ancestry-matched samples:", nrow(keep), "/", nrow(study), "\n")
fwrite(keep[, .(FID, IID)], "ancestry_keep.txt", sep=" ", col.names=FALSE)

# Covariates file: FID IID PC1..PCn
pc_cols <- paste0("PC", 1:opt$n_pcs)
covs <- keep[, c("FID","IID", pc_cols), with=FALSE]
fwrite(covs, "covariates.txt", sep="\t")
cat("Covariates written.\n")