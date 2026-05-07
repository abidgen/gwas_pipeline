#!/bin/bash
#SBATCH --job-name=verify_gwas_prs
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=128g
#SBATCH --time=48:00:00
#SBATCH --gres=lscratch:200
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

cd /home/rezaa2/Downloads/gwas-pipeline
source ~/bin/myconda
module load nextflow

export NXF_ANSI_LOG=false

nextflow run main.nf \
  -profile biowulf,conda \
  --mode full \
  --input_bed /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_all \
  --phenotype /home/rezaa2/Downloads/gwas-pipeline/data/phenotypes.csv \
  --sumstats /data/${USER}/gwas_refs/finngen_r12_lung/finngen_R12_C3_NSCLC_ADENO_EXALLC.gz \
  --local_1kg_prefix /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_all \
  --ancestry_ref_prefix /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_all \
  --assoc_test logistic \
  --pop EUR \
  --outdir results_verify_full \
  -with-report report.html \
  -with-timeline timeline.html \
  -resume