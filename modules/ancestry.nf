process ANCESTRY_PCA {
    label 'r_gwas'
    publishDir "${params.outdir}/02_ancestry", mode: 'copy'

    input:
    path(plink_files)

    output:
    path('ancestry_pass.{bed,bim,fam}'), emit: filtered_plink
    path('covariates.txt'),              emit: covariates
    path('pca/'),                        emit: pca_dir

    shell:
    '/bin/bash'

    script:
    def prefix = plink_files[0].baseName
    def ancestry_ref = params.ancestry_ref_prefix ?: params.local_1kg_prefix
    """
    mkdir -p pca

    # ── Use local 1KG reference panel ───────────────────────────────────────
    [ -n "${ancestry_ref}" ] || { echo "No ancestry reference provided. Set --ancestry_ref_prefix or --local_1kg_prefix"; exit 1; }
    cp ${ancestry_ref}.bed ref_1kg.bed
    cp ${ancestry_ref}.bim ref_1kg.bim
    cp ${ancestry_ref}.fam ref_1kg.fam
    wget -q -O 1kg.panel "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"

    [ -s ref_1kg.bed ] || { echo "ref_1kg.bed missing"; exit 1; }
    [ -s ref_1kg.bim ] || { echo "ref_1kg.bim missing"; exit 1; }
    [ -s ref_1kg.fam ] || { echo "ref_1kg.fam missing"; exit 1; }
    awk 'BEGIN{ok=1} NF!=6 {ok=0} END{exit !ok}' ref_1kg.bim || { echo "ref_1kg.bim is invalid"; head -20 ref_1kg.bim; exit 1; }

    # ── Find overlapping SNPs by CHR:BP ─────────────────────────────────────
    awk '{print \$1":"\$4}' ${prefix}.bim > study_pos.txt
    awk '{print \$1":"\$4}' ref_1kg.bim   > ref_pos.txt
    comm -12 <(sort study_pos.txt) <(sort ref_pos.txt) > overlap_pos.txt

    # ── Build extract lists from original BIMs ──────────────────────────────
    awk 'NR==FNR {keep[\$1]; next} ((\$1":"\$4) in keep) {print \$2}' overlap_pos.txt ${prefix}.bim > study_extract.txt
    awk 'NR==FNR {keep[\$1]; next} ((\$1":"\$4) in keep) {print \$2}' overlap_pos.txt ref_1kg.bim > ref_extract.txt

    plink --bfile ${prefix} --extract study_extract.txt --make-bed --out study_overlap
    plink --bfile ref_1kg   --extract ref_extract.txt   --make-bed --out ref_overlap

    # ── Force both datasets to use same variant IDs: CHR:BP ─────────────────
    awk 'BEGIN{OFS="\t"} {\$2=\$1":"\$4; print}' study_overlap.bim > study_overlap.bim.tmp && mv study_overlap.bim.tmp study_overlap.bim
    awk 'BEGIN{OFS="\t"} {\$2=\$1":"\$4; print}' ref_overlap.bim   > ref_overlap.bim.tmp   && mv ref_overlap.bim.tmp   ref_overlap.bim

    # ── Remove duplicate samples already present in study from reference ────
    awk 'NR==FNR {seen[\$1 FS \$2]=1; next} ((\$1 FS \$2) in seen) {print \$1, \$2}' study_overlap.fam ref_overlap.fam > dup_samples.txt || true

    if [ -s dup_samples.txt ]; then
        plink --bfile ref_overlap --remove dup_samples.txt --make-bed --out ref_nodup
    else
        cp ref_overlap.bed ref_nodup.bed
        cp ref_overlap.bim ref_nodup.bim
        cp ref_overlap.fam ref_nodup.fam
    fi

    # ── Merge study + reference ─────────────────────────────────────────────
    plink --bfile study_overlap --bmerge ref_nodup --make-bed --out merged || true
    if [ -f merged-merge.missnp ]; then
        plink --bfile study_overlap --exclude merged-merge.missnp --make-bed --out study_clean
        plink --bfile ref_nodup     --exclude merged-merge.missnp --make-bed --out ref_clean

        awk 'BEGIN{OFS="\t"} {\$2=\$1":"\$4; print}' study_clean.bim > study_clean.bim.tmp && mv study_clean.bim.tmp study_clean.bim
        awk 'BEGIN{OFS="\t"} {\$2=\$1":"\$4; print}' ref_clean.bim   > ref_clean.bim.tmp   && mv ref_clean.bim.tmp   ref_clean.bim

        plink --bfile study_clean --bmerge ref_clean --make-bed --out merged
    fi

    # ── LD prune + PCA ──────────────────────────────────────────────────────
    plink --bfile merged --indep-pairwise 50 5 0.2 --out pca_prune
    plink --bfile merged --extract pca_prune.prune.in --make-bed --out merged_pruned
    plink --bfile merged_pruned --pca ${params.n_pcs} --out pca/merged_pca

    # ── Ancestry filter + covariates ────────────────────────────────────────
    Rscript ${projectDir}/scripts/run_ancestry_pca.R \
        --n_pcs ${params.n_pcs} \
        --pop ${params.pop}

    # ── Subset to ancestry-matched samples ──────────────────────────────────
    plink --bfile ${prefix} --keep ancestry_keep.txt --make-bed --out ancestry_pass

    [ -s ancestry_pass.bed ] || { echo "ancestry_pass.bed missing"; exit 1; }
    [ -s ancestry_pass.bim ] || { echo "ancestry_pass.bim missing"; exit 1; }
    [ -s ancestry_pass.fam ] || { echo "ancestry_pass.fam missing"; exit 1; }
    [ -s covariates.txt ]    || { echo "covariates.txt missing"; exit 1; }
    """
}