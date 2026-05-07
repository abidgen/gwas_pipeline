#!/bin/bash
#SBATCH --job-name=dl_1kg_ref
#SBATCH --partition=norm
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g
#SBATCH --time=48:00:00
#SBATCH --gres=lscratch:300
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

cd /home/rezaa2/Downloads/gwas-pipeline
source ~/bin/myconda

# Do not run 16 PLINK imports at once with 64 GB total memory.
# 4-6 parallel chromosome jobs is a safer starting point.
export MAX_JOBS=6

# Conda env containing plink 1.9
export GWAS_PLINK_ENV=gwas_plink19

# Memory caps passed into download_1kg_ref.sh
# 6 * 8000 MB = 48 GB worst-case during parallel chr conversion
export CHR_MEM_MB=8000

# Leave headroom below 64 GB requested from Slurm
export MERGE_MEM_MB=56000

# Merge can use a few threads, but does not need all 16
export MERGE_THREADS=4

bash download_1kg_ref.sh /data/${USER}/gwas_refs/1000g_phase3_allchr