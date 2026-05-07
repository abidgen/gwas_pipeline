process DOWNLOAD_PUBLIC_DATA {
    label 'plink'
    publishDir "${params.outdir}/00_raw", mode: 'copy'

    output:
    path("1kg_all.{bed,bim,fam}"), emit: plink_files

    script:
    """
    set -euo pipefail

    # Use prebuilt local 1KG PLINK files if provided
    if [ -n "${params.local_1kg_prefix ?: ''}" ]; then
        echo "Using local 1KG reference: ${params.local_1kg_prefix}"

        [ -s "${params.local_1kg_prefix}.bed" ] || { echo "Missing ${params.local_1kg_prefix}.bed"; exit 1; }
        [ -s "${params.local_1kg_prefix}.bim" ] || { echo "Missing ${params.local_1kg_prefix}.bim"; exit 1; }
        [ -s "${params.local_1kg_prefix}.fam" ] || { echo "Missing ${params.local_1kg_prefix}.fam"; exit 1; }

        cp "${params.local_1kg_prefix}.bed" 1kg_all.bed
        cp "${params.local_1kg_prefix}.bim" 1kg_all.bim
        cp "${params.local_1kg_prefix}.fam" 1kg_all.fam
        exit 0
    fi

    mkdir -p vcf_raw

    for CHR in {1..22}; do
        wget -q -O vcf_raw/chr\${CHR}.vcf.gz \
            "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr\${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
        wget -q -O vcf_raw/chr\${CHR}.vcf.gz.tbi \
            "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr\${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi"
    done

    wget -q -O integrated_call_samples_v3.panel \
        "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"

    for CHR in {1..22}; do
        plink --vcf vcf_raw/chr\${CHR}.vcf.gz \
            --make-bed \
            --biallelic-only \
            --snps-only just-acgt \
            --vcf-half-call m \
            --out chr\${CHR}_raw
    done

    rm -f merge_list.txt
    for CHR in {2..22}; do
        echo chr\${CHR}_raw.bed chr\${CHR}_raw.bim chr\${CHR}_raw.fam >> merge_list.txt
    done

    plink --bfile chr1_raw --merge-list merge_list.txt --make-bed --out 1kg_all
    """
}