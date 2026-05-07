/*
 * MODULE: PRS
 * Calls external R scripts to avoid heredoc quoting issues.
 * scripts/run_prs.R       — detects columns, harmonises FinnGen format
 * scripts/summarise_prs.R — plots PRS distribution, computes R2
 */

process PRS {
    label 'prsice'
    publishDir "${params.outdir}/05_prs", mode: 'copy'

    input:
    path(target_plink)
    path(sumstats)
    path(ref_plink)
    path(phenotype)

    output:
    path('prs_results/'),          emit: prs_dir
    path('prs_results/*.best'),    emit: best_prs
    path('prs_results/*.summary'), emit: prs_summary

    shell:
    '/bin/bash'

    script:
    def tgt_prefix  = target_plink[0].baseName
    def ref_prefix  = ref_plink[0].baseName
    def binary_flag = params.assoc_test == 'logistic' ? '--binary-target T' : '--binary-target F'
    """
    mkdir -p prs_results

    module load prsice




    # ── Step 1: prepare phenotype (FID IID PHE) ─────────────────────────────
    awk -F',' 'NR>1 {print \$1, \$2, \$3}' ${phenotype} > prsice_pheno.txt
    [ -s prsice_pheno.txt ] || { echo "prsice_pheno.txt missing"; exit 1; }
    echo "=== prsice_pheno.txt ==="
    head -5 prsice_pheno.txt

    # ── Step 2: detect columns + harmonise sumstats via R script ────────────
    Rscript ${projectDir}/scripts/run_prs.R \
        --sumstats ${sumstats} \
        --phenotype ${phenotype} \
        --assoc_test ${params.assoc_test}

    # ── Step 3: load config written by R script ─────────────────────────────
    source sumstat_config.env

    echo "=== PRSice-2 config ==="
    echo "  Sumstats : \${SUMSTATS_FILE}"
    echo "  SNP      : \${SNP_COL}"
    echo "  CHR      : \${CHR_COL}"
    echo "  BP       : \${BP_COL}"
    echo "  A1       : \${A1_COL}"
    echo "  P        : \${P_COL}"
    echo "  Stat     : \${STAT_COL} (\${STAT_TYPE})"

    # ── Step 4: make target/reference IDs match base (CHR:BP) ──────────────
    cp ${tgt_prefix}.bed target_prs.bed
    cp ${tgt_prefix}.bim target_prs.bim
    cp ${tgt_prefix}.fam target_prs.fam
    awk 'BEGIN{OFS="\\t"} {\$2=\$1":"\$4; print}' target_prs.bim > target_prs.bim.tmp && mv target_prs.bim.tmp target_prs.bim

    cp ${ref_prefix}.bed ref_prs.bed
    cp ${ref_prefix}.bim ref_prs.bim
    cp ${ref_prefix}.fam ref_prs.fam
    awk 'BEGIN{OFS="\\t"} {\$2=\$1":"\$4; print}' ref_prs.bim > ref_prs.bim.tmp && mv ref_prs.bim.tmp ref_prs.bim

    # ── Step 5: run PRSice-2 ────────────────────────────────────────────────
    PRSice.R \
        --prsice \$(which PRSice) \
        --base   "\${SUMSTATS_FILE}" \
        --target target_prs \
        --ld     ref_prs \
        --pheno  prsice_pheno.txt \
        ${binary_flag} \
        --snp    "\${SNP_COL}" \
        --chr    "\${CHR_COL}" \
        --bp     "\${BP_COL}" \
        --A1     "\${A1_COL}" \
        --pvalue "\${P_COL}" \
        \${STAT_FLAG} \
        --bar-levels 0.001,0.01,0.05,0.1,0.2,0.3,0.4,0.5 \
        --fastscore \
        --all-score \
        --clump-r2 ${params.prs_clump_r2} \
        --clump-kb ${params.prs_clump_kb} \
        --thread ${task.cpus} \
        --out prs_results/prs

    [ -s prs_results/prs.best ]    || { echo "prs.best missing"; exit 1; }
    [ -s prs_results/prs.summary ] || { echo "prs.summary missing"; exit 1; }

    # ── Step 6: summarise + plot results ────────────────────────────────────
    Rscript ${projectDir}/scripts/summarise_prs.R \
        --best prs_results/prs.best \
        --phenotype ${phenotype} \
        --assoc_test ${params.assoc_test} \
        --stat_type "\${STAT_TYPE}" \
        --outdir prs_results
    """
}