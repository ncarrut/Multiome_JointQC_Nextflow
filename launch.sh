#!/bin/bash
#SBATCH --job-name=qc-nf
#SBATCH --time=168:00:00
#SBATCH --mem=5000M
#SBATCH --account=scjp1
#SBATCH --output=%u-%x-%j.log
#SBATCH --error=%u-%x-%j.err
#SBATCH --mail-user=xiaoouw@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --signal=B:TERM@60

module load openjdk/18.0.1.1
module load singularity/4.1.3

exec nextflow run -resume main.nf \
  --sample_id 10k_PBMC_Multiome_nextgem_Chromium_X \
  --rna_results_dir /home/xiaoouw/0_multiomePipeline/xiaoouw/results/pbmc_nf/RNA/ \
  --atac_results_dir /home/xiaoouw/0_multiomePipeline/xiaoouw/results/pbmc_nf/ATAC/ \
  --outdir /home/xiaoouw/0_multiomePipeline/xiaoouw/results/pbmc_nf/joint_QC_01092026
