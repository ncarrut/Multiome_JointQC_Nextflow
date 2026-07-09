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
module load singularity/4.4.1

exec nextflow run -resume main.nf \
  --params-file /nfs/mm-isilon/bioinfcore/ActiveProjects/ncarruth/Parker_P30/muscle_atlas/individual_datasets/Rubenstein_2025/Multiome_QC/Multiome_JointQC_Nextflow/library-config.json
  --outdir /nfs/mm-isilon/bioinfcore/ActiveProjects/ncarruth/Parker_P30/muscle_atlas/individual_datasets/Rubenstein_2025/Multiome_QC/results
