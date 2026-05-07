process ASSOCIATION {
    label 'r_gwas'
    publishDir "${params.outdir}/03_association", mode: 'copy'

    input:
    path(plink_files)
    path(phenotype)
    path(covariates)

    output:
    path('gwas_results.assoc'), emit: sumstats
    path('gwas_results.*'),     emit: all_results

    script:
    def prefix     = plink_files[0].baseName
    def test_flag  = params.assoc_test == 'logistic' ? '--logistic hide-covar' : '--linear hide-covar'
    def assoc_file = params.assoc_test == 'logistic' ? 'gwas_raw.assoc.logistic' : 'gwas_raw.assoc.linear'
    """
    # ── Prepare phenotype + covariates ────────────────────────────────────────
    Rscript ${projectDir}/scripts/run_association.R \
        ${phenotype} \
        ${covariates} \
        ${params.gwas_pval}

    [ -s plink_pheno.txt ]    || { echo "plink_pheno.txt missing"; exit 1; }
    [ -s all_covariates.txt ] || { echo "all_covariates.txt missing"; exit 1; }

    echo "=== plink_pheno.txt ==="
    head -5 plink_pheno.txt
    echo "=== all_covariates.txt ==="
    head -5 all_covariates.txt

    # ── Run association test ──────────────────────────────────────────────────
    plink \
        --bfile ${prefix} \
        --pheno plink_pheno.txt \
        --covar all_covariates.txt \
        --allow-no-sex \
        ${test_flag} \
        --ci 0.95 \
        --pfilter 1 \
        --out gwas_raw

    [ -s ${assoc_file} ] || { echo "${assoc_file} missing"; exit 1; }
    cp ${assoc_file} gwas_results.assoc

    # ── Lambda + significant hits ─────────────────────────────────────────────
    Rscript ${projectDir}/scripts/summarise_association.R \
        gwas_results.assoc \
        ${params.gwas_pval}
    """
}