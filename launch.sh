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
module load singularity/4.1.5

exec nextflow run -resume Multiome_JointQC_Nextflow/main.nf \
  --samplesheet /nfs/mm-isilon/bioinfcore/ActiveProjects/Lukacs_nlukacs_CU2-Multiome_ncarruth_rcavalca_12540-AE/analysis_from_cell_ranger/Multiome_JointQC_Nextflow/samplesheet.txt
  