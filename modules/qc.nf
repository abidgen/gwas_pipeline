/*
 * MODULE: QC
 * R code moved to scripts/run_qc.R to avoid Nextflow variable conflicts.
 */

process QC {
    label 'r_gwas'
    publishDir "${params.outdir}/01_qc", mode: 'copy'

    input:
    path(plink_files)

    output:
    path('qc_pass.{bed,bim,fam}'), emit: clean_plink
    path('qc_stats/'),             emit: qc_stats

    script:
    def prefix = plink_files[0].baseName
    """
    mkdir -p qc_stats

    # ── 1. Initial missingness ────────────────────────────────────────────────
    plink --bfile ${prefix} --missing --out qc_stats/initial_missing

    # ── 2. Sample call rate ───────────────────────────────────────────────────
    plink --bfile ${prefix} --mind ${params.mind} --make-bed --out step1_mind

    # ── 3. Sex check (skipped if no X chromosome SNPs present) ───────────────
    N_XCHR=\$(awk '\$1==23' step1_mind.bim | wc -l)
    if [ "\$N_XCHR" -gt 0 ]; then
        plink --bfile step1_mind --check-sex --out qc_stats/sexcheck
        awk 'NR>1 && \$5=="PROBLEM" {print \$1, \$2}' qc_stats/sexcheck.sexcheck \
            > sex_fails.txt || touch sex_fails.txt
        echo "Sex check run on \$N_XCHR X-chr SNPs. Failures: \$(wc -l < sex_fails.txt)"
    else
        echo "No X chromosome SNPs found - skipping sex check (normal for chr22-only demo)" \
            > qc_stats/sexcheck_skipped.txt
        touch sex_fails.txt
    fi
    plink --bfile step1_mind --remove sex_fails.txt --make-bed --out step2_sex

    # ── 4. Heterozygosity (via external R script) ─────────────────────────────
    plink --bfile step2_sex --het --out qc_stats/het
    Rscript ${projectDir}/scripts/run_qc.R
    plink --bfile step2_sex --remove het_fails.txt --make-bed --out step3_het

    # ── 5. Variant QC ─────────────────────────────────────────────────────────
    plink --bfile step3_het \
        --geno ${params.geno} \
        --maf  ${params.maf}  \
        --hwe  ${params.hwe}  \
        --make-bed --out step4_snpqc

    # ── 6. Relatedness ────────────────────────────────────────────────────────
    plink --bfile step4_snpqc --indep-pairwise 50 5 0.2 --out ld_prune_list
    plink --bfile step4_snpqc \
        --extract ld_prune_list.prune.in \
        --genome --min ${params.pihat} \
        --out qc_stats/ibd
    awk 'NR>1 {print \$1, \$2}' qc_stats/ibd.genome \
        | sort -u > related_remove.txt || touch related_remove.txt
    plink --bfile step4_snpqc --remove related_remove.txt --make-bed --out qc_pass

    # ── 7. Summary ────────────────────────────────────────────────────────────
    plink --bfile qc_pass --missing --out qc_stats/final_missing
    echo "SNPs remaining: \$(wc -l < qc_pass.bim)"     > qc_stats/qc_summary.txt
    echo "Samples remaining: \$(wc -l < qc_pass.fam)" >> qc_stats/qc_summary.txt
    cat qc_stats/qc_summary.txt
    """
}