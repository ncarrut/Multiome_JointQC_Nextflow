# Multiome Joint QC (Nextflow)

This workflow runs the `joint_qc.py` script to generate joint RNA/ATAC QC plots.

## Requirements

- Nextflow (DSL2)
- Singularity/Apptainer module available on compute nodes
- Container image (default: `docker://ncarrut/singlecell_qc:second`)

## Inputs

Required parameters:
- `--params-file` : json file listing sample RNA and ATAC data locations for each sample (see library-config.json)
- `--outdir` : output directory for results (default: `results`)
- `--filter_MT_ATAC` : whether to filter ATAC nuclei on a %chrMT threshold (default: `false`)

## Run

```bash
nextflow run main.nf \
  --params-file library-config.json \
  --outdir /path/to/output
```

## Outputs

Results are written to `--outdir/`:
- `<sample_id>_qcPlot.png`
- `<sample_id>_upsetPlot.png`
- `<sample_id>_metrics.txt`
- `<sample_id>_joint_qc.log`

## Container

The container is configured in `nextflow.config`:

```
process.container = 'library://alicewang24/python/jointqc_20260113:latest'
```

If you want to use a local `.sif`, replace the value with the absolute path and
ensure the directory is bound in `containerOptions`.
