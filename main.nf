#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.mode         = 'full'
params.input_bed    = null
params.sumstats     = null
params.phenotype    = 'data/phenotypes.csv'
params.outdir       = 'results'
params.ref_panel    = null
params.gwas_pval    = 5e-8
params.prs_pval     = 0.5
params.pop          = 'EUR'

include { DOWNLOAD_PUBLIC_DATA } from './modules/download'
include { QC                   } from './modules/qc'
include { ANCESTRY_PCA         } from './modules/ancestry'
include { ASSOCIATION          } from './modules/association'
include { PLOTS                } from './modules/plots'
include { PRS                  } from './modules/prs'

workflow {
    log.info """
    ╔══════════════════════════════════════════╗
    ║         GWAS NEXTFLOW PIPELINE           ║
    ║  mode      : ${params.mode}           
    ║  outdir    : ${params.outdir}
    ║  phenotype : ${params.phenotype}
    ║  gwas_pval : ${params.gwas_pval}
    ╚══════════════════════════════════════════╝
    """.stripIndent()

    pheno_ch = Channel.fromPath(params.phenotype, checkIfExists: true)

    if (params.input_bed) {
        bed_ch = Channel
            .fromFilePairs("${params.input_bed}.{bed,bim,fam}", size: 3)
            .map { prefix, files -> files }
    }
    else {
        DOWNLOAD_PUBLIC_DATA()
        bed_ch = DOWNLOAD_PUBLIC_DATA.out.plink_files
    }

    QC(bed_ch)
    ANCESTRY_PCA(QC.out.clean_plink)

    if (params.mode in ['full', 'assoc']) {
        ASSOCIATION(
            ANCESTRY_PCA.out.filtered_plink,
            pheno_ch,
            ANCESTRY_PCA.out.covariates
        )
        PLOTS(ASSOCIATION.out.sumstats)
    }

    if (params.mode in ['full', 'prs']) {
        def sumstats_ch = params.sumstats
            ? Channel.fromPath(params.sumstats, checkIfExists: true)
            : ASSOCIATION.out.sumstats

        def prs_ref_prefix = params.ref_panel ?: params.local_1kg_prefix
        assert prs_ref_prefix : "Provide --ref_panel or --local_1kg_prefix for PRS reference"

        def ref_ch = Channel
            .fromFilePairs("${prs_ref_prefix}.{bed,bim,fam}", size: 3)
            .map { p, f -> f }

        PRS(
            ANCESTRY_PCA.out.filtered_plink,
            sumstats_ch,
            ref_ch,
            pheno_ch
        )
    }
}

workflow.onComplete {
    log.info "Pipeline complete! Results in: ${params.outdir}"
}