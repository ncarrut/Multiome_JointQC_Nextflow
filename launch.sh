#!/bin/bash
#SBATCH --job-name=qc-nf
#SBATCH --time=168:00:00
#SBATCH --mem=5000M
#SBATCH --account=<account>
#SBATCH --output=%u-%x-%j.log
#SBATCH --error=%u-%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --signal=B:TERM@60

module load openjdk/18.0.1.1
module load singularity/4.4.1

BASE_DIR=/nfs/mm-isilon/bioinfcore/ActiveProjects/ncarruth/Parker_P30/muscle_atlas/individual_datasets/Kedlian_2024/Multiome_QC/Multiome_JointQC_Nextflow

nextflow run "${BASE_DIR}/main.nf" \
    --samplesheet "${BASE_DIR}/library-config.tsv" \
    --outdir "${BASE_DIR}/../results"

