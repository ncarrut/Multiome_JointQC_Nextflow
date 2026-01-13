# Multiome Joint QC (Nextflow)

This workflow runs the `joint_qc.py` script to generate joint RNA/ATAC QC plots
and metrics for a single sample.

## Requirements

- Nextflow (DSL2)
- Singularity/Apptainer module available on compute nodes
- Container image (default: `library://alicewang24/python/jointqc_20260113:latest`)

## Inputs

Required parameters:
- `--sample_id` : sample name
- `--rna_results_dir` : RNA results directory (output of https://github.com/porchard/snRNAseq-NextFlow)
- `--atac_results_dir` : ATAC results directory (output of https://github.com/porchard/snATACseq-NextFlow)
- `--outdir` : output directory for results (default: `results`)

## Run

```bash
nextflow run main.nf \
  --sample_id 10k_PBMC_Multiome_nextgem_Chromium_X \
  --rna_results_dir /path/to/RNA/ \
  --atac_results_dir /path/to/ATAC/ \
  --outdir /path/to/output
```

## Outputs

Results are written to `--outdir/<sample_id>/`:
- `qcPlot.png`
- `upsetPlot.png`
- `metrics.txt`

## Container

The container is configured in `nextflow.config`:

```
process.container = 'library://alicewang24/python/jointqc_20260113:latest'
```

If you want to use a local `.sif`, replace the value with the absolute path and
ensure the directory is bound in `containerOptions`.
