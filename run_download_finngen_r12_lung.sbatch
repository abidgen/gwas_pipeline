#!/bin/bash
#SBATCH --job-name=dl_finngen_r12
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

export MAX_JOBS=${SLURM_CPUS_PER_TASK}

bash download_finngen_r12_lung.sh /data/${USER}/gwas_refs/finngen_r12_lung