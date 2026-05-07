process PLOTS {
    label 'r_gwas'
    publishDir "${params.outdir}/04_plots", mode: 'copy'
    input:
    path(sumstats)
    output:
    path("manhattan.png"), emit: manhattan
    path("qqplot.png"),    emit: qqplot
    path("lambda.txt"),    emit: lambda
    script:
    """
    Rscript ${projectDir}/scripts/plot_manhattan_qq.R \
        --sumstats ${sumstats} --pval_col P --chr_col CHR --pos_col BP --snp_col SNP \
        --gwas_sig ${params.gwas_pval} --out_manhattan manhattan.png \
        --out_qq qqplot.png --out_lambda lambda.txt
    """
}
