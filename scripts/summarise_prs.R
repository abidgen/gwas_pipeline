#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(optparse)
})

opt_list <- list(
    make_option("--best",       type = "character"),
    make_option("--phenotype",  type = "character"),
    make_option("--assoc_test", type = "character", default = "logistic"),
    make_option("--stat_type",  type = "character", default = "OR"),
    make_option("--outdir",     type = "character", default = "prs_results")
)

parser <- OptionParser(option_list = opt_list)
opt <- parse_args(parser)

if (is.null(opt$best) || is.null(opt$phenotype)) {
    print_help(parser)
    stop("Both --best and --phenotype are required.")
}

best  <- fread(opt$best)
pheno <- fread(opt$phenotype)

df <- merge(best, pheno[, .(FID, IID, PHE)], by = c("FID", "IID"))

# Recode PLINK phenotype: 1=control, 2=case -> 0/1 for glm
df[, PHE01 := fifelse(PHE == 2, 1,
               fifelse(PHE == 1, 0, NA_real_))]
df <- df[!is.na(PHE01)]

df[, PHE_LABEL := factor(PHE, levels = c(1, 2), labels = c("Control", "Case"))]

if (opt$assoc_test == "logistic") {
    fit0   <- glm(PHE01 ~ 1,   data = df, family = binomial())
    fit1   <- glm(PHE01 ~ PRS, data = df, family = binomial())
    r2     <- as.numeric(1 - logLik(fit1) / logLik(fit0))
    metric <- sprintf("Nagelkerke R2 = %.4f", r2)
} else {
    fit    <- lm(PHE ~ PRS, data = df)
    r2     <- summary(fit)$r.squared
    metric <- sprintf("R2 = %.4f", r2)
}

cat(metric, "\n")
cat(sprintf("Effect size format: %s\n", opt$stat_type))

p <- ggplot(df, aes(PRS, fill = PHE_LABEL)) +
    geom_density(alpha = 0.6) +
    labs(
        title    = paste("PRS distribution -", metric),
        subtitle = paste("Effect size format:", opt$stat_type),
        x        = "Polygenic Risk Score",
        fill     = "Phenotype"
    ) +
    theme_bw()

ggsave(file.path(opt$outdir, "prs_distribution.png"), p, width = 8, height = 5, dpi = 150)

df[, quartile := cut(PRS, quantile(PRS, 0:4/4), include.lowest = TRUE, labels = 1:4)]

qs <- df[, .(
    mean_prs = mean(PRS),
    n = .N,
    cases = sum(PHE01 == 1),
    controls = sum(PHE01 == 0)
), by = quartile]

fwrite(qs, file.path(opt$outdir, "prs_quartile_summary.txt"), sep = "\t")

cat("PRS summary complete.\n")