#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(qqman); library(optparse) })
opt_list <- list(
    make_option("--sumstats",      type="character"),
    make_option("--pval_col",      type="character", default="P"),
    make_option("--chr_col",       type="character", default="CHR"),
    make_option("--pos_col",       type="character", default="BP"),
    make_option("--snp_col",       type="character", default="SNP"),
    make_option("--gwas_sig",      type="double",    default=5e-8),
    make_option("--out_manhattan", type="character", default="manhattan.png"),
    make_option("--out_qq",        type="character", default="qqplot.png"),
    make_option("--out_lambda",    type="character", default="lambda.txt")
)
opt <- parse_args(OptionParser(option_list=opt_list))
res <- fread(opt$sumstats)
setnames(res, c(opt$chr_col,opt$pos_col,opt$snp_col,opt$pval_col), c("CHR","BP","SNP","P"))
res <- res[!is.na(P) & P>0 & CHR %in% 1:22]
res[, CHR := as.integer(CHR)]
chisq  <- qchisq(res$P, df=1, lower.tail=FALSE)
lambda <- median(chisq) / qchisq(0.5, df=1, lower.tail=FALSE)
writeLines(sprintf("lambda=%.4f", lambda), opt$out_lambda)
png(opt$out_manhattan, width=1800, height=700, res=150)
manhattan(res, chr="CHR", bp="BP", snp="SNP", p="P",
          genomewideline=-log10(opt$gwas_sig), suggestiveline=-log10(1e-5),
          main=sprintf("Manhattan plot (λ=%.3f)", lambda), col=c("#3B6D11","#0C447C"))
dev.off()
png(opt$out_qq, width=700, height=700, res=150)
qq(res$P, main=sprintf("QQ plot (λ=%.3f)", lambda))
dev.off()
cat("Plots written.\n")
