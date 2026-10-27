#!/bin/bash
#SBATCH --job-name=qc-nf
#SBATCH --time=168:00:00
#SBATCH --mem=5000M
#SBATCH --account=scjp1
#SBATCH --output=%u-%x-%j.log
#SBATCH --error=%u-%x-%j.err
#SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=uniqname@umich.edu
#SBATCH --signal=B:TERM@60

# Usage (from the directory where results/ and work/ should be written):
#   sbatch /path/to/Multiome_JointQC_Nextflow/conf/greatlakes/launch.sh samplesheet.txt [results_dir]

set -euo pipefail

SAMPLESHEET="${1:?usage: launch.sh <samplesheet.txt> [results_dir]}"
RESULTS="${2:-results}"

# Locate the pipeline checkout. Under sbatch the script runs from a spool copy,
# so ask slurm for the original script path. Set PIPELINE_DIR to override.
if [[ -z "${PIPELINE_DIR:-}" ]]; then
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
        SCRIPT_PATH="$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/ Command=/{print $2}' | awk '{print $1}')"
    else
        SCRIPT_PATH="${BASH_SOURCE[0]}"
    fi
    PIPELINE_DIR="$(cd "$(dirname "$SCRIPT_PATH")/../.." && pwd)"
fi

module load openjdk/18.0.1.1
module load singularity/4.4.1
module load Bioinformatics
module load nextflow

exec nextflow run -resume "${PIPELINE_DIR}/main.nf" \
  -profile greatlakes \
  --samplesheet "${SAMPLESHEET}" \
  --results "${RESULTS}"
