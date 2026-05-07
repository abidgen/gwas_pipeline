#!/usr/bin/env Rscript
# =============================================================================
# validate_simulation.R
# Checks whether the pipeline recovered the known causal SNPs from simulation.
# Produces a validation report + plots.
#
# Usage:
#   Rscript scripts/validate_simulation.R \
#       --gwas_results  results/03_association/gwas_results.assoc \
#       --prs_best      results/05_prs/prs_results/prs.best \
#       --causal_snps   simulated_ld_data/causal_snps.txt \
#       --phenotype     data/phenotypes.csv \
#       --outdir        results/validation/
# =============================================================================
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(optparse)
})

opt_list <- list(
    make_option("--gwas_results", type="character", help="GWAS assoc file"),
    make_option("--prs_best",     type="character", help="PRSice .best file"),
    make_option("--causal_snps",  type="character", help="True causal SNP list"),
    make_option("--phenotype",    type="character", help="Phenotype CSV"),
    make_option("--outdir",       type="character", default="results/validation")
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

cat("=======================================================\n")
cat("  GWAS Pipeline Validation Report\n")
cat("=======================================================\n\n")

# ── 1. Load data ──────────────────────────────────────────────────────────────
causal <- fread(opt$causal_snps, header=FALSE, col.names="SNP")
cat(sprintf("True causal SNPs: %d\n\n", nrow(causal)))

# ── 2. GWAS recovery check ────────────────────────────────────────────────────
if (!is.null(opt$gwas_results) && file.exists(opt$gwas_results)) {
    cat("--- GWAS Association Recovery ---\n")
    gwas <- fread(opt$gwas_results)
    gwas <- gwas[!is.na(P)]

    # Lambda
    chisq  <- qchisq(gwas$P, df=1, lower.tail=FALSE)
    lambda <- median(chisq) / qchisq(0.5, df=1, lower.tail=FALSE)
    cat(sprintf("Genomic inflation (lambda): %.4f", lambda))
    if (lambda < 1.05) cat("  ✓ GOOD (< 1.05)\n")
    else if (lambda < 1.10) cat("  ~ ACCEPTABLE (1.05-1.10)\n")
    else cat("  ✗ INFLATED (> 1.10) — check population stratification\n")

    # How many causal SNPs are in the GWAS results
    causal_in_gwas <- gwas[SNP %in% causal$SNP]
    cat(sprintf("\nCausal SNPs present in GWAS results: %d / %d\n",
        nrow(causal_in_gwas), nrow(causal)))

    # Rank of causal SNPs
    gwas_ranked <- gwas[order(P)]
    gwas_ranked[, rank := .I]
    causal_ranks <- gwas_ranked[SNP %in% causal$SNP, .(SNP, P, rank, .N)]
    if (nrow(causal_ranks) > 0) {
        cat("\nCausal SNP ranks in GWAS (lower = better recovered):\n")
        print(causal_ranks[order(rank)], row.names=FALSE)
    }

    # Recovery at different p-value thresholds
    thresholds <- c(5e-8, 1e-6, 1e-5, 1e-4, 0.05)
    cat("\nCausal SNP recovery by p-value threshold:\n")
    cat(sprintf("  %-12s  %-10s  %-10s  %-8s\n",
        "Threshold", "Recovered", "Total hits", "Power"))
    for (thr in thresholds) {
        sig_snps  <- gwas[P < thr, SNP]
        recovered <- sum(causal$SNP %in% sig_snps)
        total_sig <- length(sig_snps)
        power     <- recovered / nrow(causal)
        cat(sprintf("  %-12s  %-10d  %-10d  %.1f%%\n",
            formatC(thr, format="e", digits=0),
            recovered, total_sig, power*100))
    }

    # Manhattan plot highlighting causal SNPs
    cat("\nGenerating Manhattan plot with causal SNPs highlighted...\n")
    gwas[, is_causal := SNP %in% causal$SNP]
    gwas[, CHR := as.integer(CHR)]
    gwas <- gwas[CHR %in% 1:22 & P > 0]

    # Chromosome offsets for plotting
    chr_sizes <- gwas[, .(max_bp=max(BP)), by=CHR][order(CHR)]
    chr_sizes[, offset := cumsum(shift(max_bp, fill=0))]
    gwas <- merge(gwas, chr_sizes[,.(CHR, offset)], by="CHR")
    gwas[, pos_global := BP + offset]

    p_manh <- ggplot(gwas, aes(pos_global, -log10(P))) +
        geom_point(aes(color=factor(CHR %% 2)), size=0.4, alpha=0.5, show.legend=FALSE) +
        geom_point(data=gwas[is_causal==TRUE],
                   aes(pos_global, -log10(P)),
                   color="red", size=3, shape=18) +
        geom_hline(yintercept=-log10(5e-8), linetype="dashed", color="red",    linewidth=0.5) +
        geom_hline(yintercept=-log10(1e-5),  linetype="dashed", color="orange", linewidth=0.5) +
        scale_color_manual(values=c("#3B6D11","#0C447C")) +
        scale_x_continuous(
            breaks = chr_sizes[, offset + max_bp/2],
            labels = chr_sizes$CHR
        ) +
        labs(
            title    = sprintf("Manhattan plot — causal SNPs highlighted (lambda=%.3f)", lambda),
            subtitle = sprintf("Red diamonds = true causal SNPs (n=%d)", nrow(causal)),
            x = "Chromosome", y = expression(-log[10](p))
        ) +
        theme_bw(base_size=11) +
        theme(axis.text.x=element_text(size=7), panel.grid=element_blank())

    ggsave(file.path(opt$outdir, "manhattan_validated.png"),
           p_manh, width=14, height=5, dpi=150)

    # P-value distribution of causal vs non-causal SNPs
    p_enrich <- ggplot() +
        geom_histogram(data=gwas[is_causal==FALSE],
                       aes(P, fill="Non-causal"),
                       bins=50, alpha=0.6) +
        geom_histogram(data=gwas[is_causal==TRUE],
                       aes(P, fill="Causal"),
                       bins=20, alpha=0.8) +
        scale_fill_manual(values=c("Causal"="red","Non-causal"="steelblue")) +
        labs(title="P-value distribution: causal vs non-causal SNPs",
             x="P-value", y="Count", fill="") +
        theme_bw()
    ggsave(file.path(opt$outdir, "pvalue_enrichment.png"),
           p_enrich, width=8, height=5, dpi=150)

    cat(sprintf("Plots saved to %s/\n", opt$outdir))
}

# ── 3. PRS validation ─────────────────────────────────────────────────────────
if (!is.null(opt$prs_best) && file.exists(opt$prs_best)) {
    cat("\n--- PRS Validation ---\n")
    prs   <- fread(opt$prs_best)
    pheno <- fread(opt$phenotype)
    df    <- merge(prs, pheno[,.(FID,IID,PHE)], by=c("FID","IID"))

    # AUC
    df[, PHE_bin := as.integer(PHE == 2)]
    fit   <- glm(PHE_bin ~ PRS, data=df, family=binomial)
    probs <- predict(fit, type="response")

    # Simple AUC via rank sum
    pos <- probs[df$PHE_bin==1]
    neg <- probs[df$PHE_bin==0]
    auc <- mean(outer(pos, neg, ">")) + 0.5*mean(outer(pos, neg, "=="))
    cat(sprintf("PRS AUC           : %.4f", auc))
    if (auc > 0.6)       cat("  ✓ GOOD signal\n")
    else if (auc > 0.55) cat("  ~ WEAK signal\n")
    else                 cat("  ✗ No signal (AUC ~ 0.5 = random)\n")

    # Nagelkerke R2
    fit0 <- glm(PHE_bin ~ 1, data=df, family=binomial)
    r2   <- as.numeric(1 - logLik(fit)/logLik(fit0))
    cat(sprintf("Nagelkerke R2     : %.4f\n", r2))

    # Odds ratio top vs bottom quartile
    df[, quartile := cut(PRS, quantile(PRS, 0:4/4),
                         include.lowest=TRUE, labels=1:4)]
    q1 <- df[quartile==1]; q4 <- df[quartile==4]
    or_data <- rbind(q1, q4)
    fit_or  <- glm(PHE_bin ~ quartile, data=or_data, family=binomial)
    or_val  <- exp(coef(fit_or)[2])
    cat(sprintf("OR (Q4 vs Q1)     : %.3f\n", or_val))
    cat(sprintf("Expected OR range : 1.5-3.0 for h2=%.2f, N=%d\n",
        0.1, nrow(df)))

    # PRS distribution plot by case/control
    p_prs <- ggplot(df, aes(PRS, fill=factor(PHE))) +
        geom_density(alpha=0.6) +
        scale_fill_manual(values=c("1"="#185FA5","2"="#E24B4A"),
                          labels=c("1"="Control","2"="Case")) +
        labs(
            title    = sprintf("PRS distribution — AUC=%.3f, R2=%.4f", auc, r2),
            subtitle = "Case curve should be shifted RIGHT of control curve",
            x = "Polygenic Risk Score", fill = ""
        ) +
        theme_bw()
    ggsave(file.path(opt$outdir, "prs_validated.png"),
           p_prs, width=8, height=5, dpi=150)

    # Quartile bar plot
    q_summary <- df[, .(
        cases    = sum(PHE==2),
        controls = sum(PHE==1),
        pct_case = mean(PHE==2)*100
    ), by=quartile][order(quartile)]

    p_qbar <- ggplot(q_summary, aes(quartile, pct_case, fill=quartile)) +
        geom_col(show.legend=FALSE) +
        geom_text(aes(label=sprintf("%.1f%%", pct_case)), vjust=-0.5, size=3.5) +
        scale_fill_manual(values=c("#C6DBF0","#7EB4DA","#2171B5","#08306B")) +
        labs(
            title    = "Case proportion by PRS quartile",
            subtitle = "Should increase from Q1 to Q4 for a predictive PRS",
            x = "PRS Quartile", y = "% Cases"
        ) +
        ylim(0, max(q_summary$pct_case)*1.2) +
        theme_bw()
    ggsave(file.path(opt$outdir, "prs_quartiles_validated.png"),
           p_qbar, width=6, height=5, dpi=150)
}

# ── 4. Final summary ──────────────────────────────────────────────────────────
cat("\n=======================================================\n")
cat("  Validation Output Files\n")
cat("=======================================================\n")
files <- list.files(opt$outdir, full.names=FALSE)
for (f in files) cat(sprintf("  %s/%s\n", opt$outdir, f))
cat("\nValidation complete.\n")