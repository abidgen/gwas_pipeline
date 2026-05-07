#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(optparse)
})

opt_list <- list(
    make_option("--sumstats",   type="character"),
    make_option("--phenotype",  type="character"),
    make_option("--assoc_test", type="character", default="logistic"),
    make_option("--outdir",     type="character", default="prs_results")
)
opt <- parse_args(OptionParser(option_list=opt_list))

# ── Step 1: detect columns and effect size type ───────────────────────────────
if (grepl("gz$", opt$sumstats)) {
    hdr <- fread(cmd = paste("zcat", shQuote(opt$sumstats)), nrows = 5)
} else {
    hdr <- fread(opt$sumstats, nrows = 5)
}

cols <- tolower(names(hdr))
cat("Detected columns:", paste(names(hdr), collapse = ", "), "\n")

snp_col <- names(hdr)[which(cols %in% c("snp","rsid","rsids","variant","snpid","markername"))[1]]
if (is.na(snp_col)) snp_col <- names(hdr)[grep("rs|snp|variant", cols)[1]]
chr_col <- names(hdr)[which(cols %in% c("chr","chrom","chromosome","#chrom"))[1]]
bp_col  <- names(hdr)[which(cols %in% c("bp","pos","position","chromstart"))[1]]
a1_col  <- names(hdr)[which(cols %in% c("a1","alt","effect_allele","ea","allele1"))[1]]
p_col   <- names(hdr)[which(cols %in% c("p","pval","p_value","p.value","pvalue"))[1]]

has_or   <- any(cols %in% c("or","odds_ratio"))
has_beta <- any(cols %in% c("beta","b","effect"))

if (has_or && !has_beta) {
    stat_col  <- names(hdr)[which(cols %in% c("or","odds_ratio"))[1]]
    stat_type <- "OR"
    stat_flag <- "--or"
} else if (has_beta) {
    stat_col  <- names(hdr)[which(cols %in% c("beta","b","effect"))[1]]
    stat_type <- "beta"
    stat_flag <- "--beta"
} else {
    stop("Cannot detect effect size column. Columns: ", paste(names(hdr), collapse = ", "))
}

cat(sprintf("Effect size: %s (%s)\n", stat_col, stat_type))
cat(sprintf("SNP=%s CHR=%s BP=%s A1=%s P=%s\n", snp_col, chr_col, bp_col, a1_col, p_col))

# ── Step 2: FinnGen harmonisation / standardisation ───────────────────────────
is_finngen <- grepl("finngen", tolower(opt$sumstats)) ||
              any(grepl("^#chrom$", tolower(names(hdr))))

if (grepl("gz$", opt$sumstats)) {
    dat <- fread(cmd = paste("zcat", shQuote(opt$sumstats)))
} else {
    dat <- fread(opt$sumstats)
}

setnames(dat, names(dat), tolower(gsub("^#", "", names(dat))))

if (is_finngen) {
    cat("FinnGen format detected - harmonising...\n")
}

# standardize key columns
if ("rsids" %in% names(dat)) setnames(dat, "rsids", "rsids_orig")

chr_col <- "chrom"
bp_col  <- "pos"
a1_col  <- "alt"
p_col   <- "pval"

# clean chromosome first
dat[, chrom := as.integer(sub("^chr", "", as.character(chrom)))]
dat <- dat[!is.na(chrom) & chrom %in% 1:22]

# force base SNP IDs to match PLINK IDs we will set in target/reference
dat[, SNP := paste0(chrom, ":", pos)]
snp_col <- "SNP"

# clean SNP ids
dat[, SNP := trimws(as.character(SNP))]
dat <- dat[!is.na(SNP) & SNP != "" & SNP != "."]

# numeric p-values
dat[, pval := suppressWarnings(as.numeric(pval))]
dat <- dat[!is.na(pval)]

# keep the most significant row per SNP
setorderv(dat, c("SNP", "pval"), c(1, 1))
before_n <- nrow(dat)
dat <- dat[, .SD[1], by = SNP]
after_n <- nrow(dat)

cat(sprintf("Deduplicated SNP IDs: %d -> %d rows\n", before_n, after_n))

fwrite(dat, "sumstats_harmonised.tsv.gz", sep = "\t", compress = "gzip")
sumstats_file <- "sumstats_harmonised.tsv.gz"
cat("Harmonised file written.\n")

# ── Step 3: write config for shell ───────────────────────────────────────────
cfg <- c(
    paste0('SUMSTATS_FILE="', sumstats_file, '"'),
    paste0('SNP_COL="',  snp_col, '"'),
    paste0('CHR_COL="',  chr_col, '"'),
    paste0('BP_COL="',   bp_col, '"'),
    paste0('A1_COL="',   a1_col, '"'),
    paste0('P_COL="',    p_col, '"'),
    paste0('STAT_COL="', stat_col, '"'),
    paste0('STAT_TYPE="', stat_type, '"'),
    paste0('STAT_FLAG="', stat_flag, '"')
)
writeLines(cfg, "sumstat_config.env")
cat("Config written to sumstat_config.env\n")