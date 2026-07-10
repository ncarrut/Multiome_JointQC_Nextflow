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

BASE_DIR=/nfs/mm-isilon/bioinfcore/ActiveProjects/ncarruth/Parker_P30/muscle_atlas/individual_datasets/Kedlian_2024/Multiome_QC/Multiome_JointQC_Nextflow

# Sample IDs repeat across data types/chemistries (e.g. donor 362C was profiled
# as both scRNA and snRNA), so each config is run separately into its own
# outdir to avoid output filename collisions.
for config in single_cell_v2 single_cell_v3 single_nuclei_v2 single_nuclei_v3 single_nuclei_mixed; do
    nextflow run "${BASE_DIR}/main.nf" \
        -params-file "${BASE_DIR}/library-config_${config}.json" \
        --outdir "${BASE_DIR}/../results/${config}"
done
