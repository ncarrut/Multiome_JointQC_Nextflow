# Multiome / scRNA / snRNA Joint QC (Nextflow)

This workflow runs QC for three assay types, selected per-sample via an `"assay"`
field in the library config:

- `multiome` (default if `assay` is omitted) — joint RNA/ATAC QC via `joint_qc.py`.
  Requires both `RNA_results_dir` and `ATAC_results_dir`.
- `snRNA` — RNA-only QC via `rna_qc.py`, using the standard exon/full-gene-body
  intron-retention ratio filter (same logic as multiome). Requires `RNA_results_dir`.
- `scRNA` — RNA-only QC via `rna_qc.py`, using **HELM**
  (`rna_fraction_mitochondrial * (1 - rna_exon_to_full_gene_body_ratio)`) in place
  of the exon-ratio filter, since whole-cell preps don't carry the same nuclear
  pre-mRNA contamination signature that the exon-ratio filter targets. Requires
  `RNA_results_dir`.

A single library-config.json can mix all three assay types in one run — each
sample is routed to the matching process automatically. However, when the same
donor was profiled more than once (e.g. as both scRNA and snRNA, or across
different chemistry versions), each of those runs needs its own config file and
`--outdir`, since the sample ID repeats and would otherwise overwrite outputs
from the other run. See `library-config_single_cell_v2.json`,
`library-config_single_cell_v3.json`, `library-config_single_nuclei_v2.json`,
`library-config_single_nuclei_v3.json`, and `library-config_single_nuclei_mixed.json`
for this project's actual per-type/chemistry configs, and `launch.sh` for how
they're each run into a separate `results/<config>/` subdirectory in one job.

## Requirements

- Nextflow (DSL2)
- Singularity/Apptainer module available on compute nodes
- Container image (default: `docker://ncarrut/singlecell_qc:second`)

## Inputs

Required parameters:
- `--params-file` : json file listing sample RNA (and, for multiome, ATAC) data locations for each sample (see `library-config_single_cell_v2.json` for an example)
- `--outdir` : output directory for results (default: `results`)
- `--filter_MT_ATAC` : whether to filter ATAC nuclei on a %chrMT threshold; multiome only (default: `false`)

Each entry under `"libraries"` in the params file takes:
- `"assay"` : one of `multiome`, `snRNA`, `scRNA`. Optional, defaults to `multiome`.
- `"RNA_results_dir"` : path to the upstream snRNAseq-NextFlow results directory. Always required.
- `"ATAC_results_dir"` : path to the upstream snATACseq-NextFlow results directory. Required only when `assay` is `multiome`.

```json
{
  "libraries": {
    "Ctrl-1_Pre": {
      "assay": "multiome",
      "RNA_results_dir": "/path/to/snRNAseq-NextFlow/results/",
      "ATAC_results_dir": "/path/to/snATACseq-NextFlow/results/"
    },
    "SampleB": {
      "assay": "snRNA",
      "RNA_results_dir": "/path/to/snRNAseq-NextFlow/results/"
    },
    "SampleC": {
      "assay": "scRNA",
      "RNA_results_dir": "/path/to/snRNAseq-NextFlow/results/"
    }
  }
}
```

## Run

```bash
nextflow run main.nf \
  --params-file library-config_single_cell_v2.json \
  --outdir results/single_cell_v2
```

To run every config for this project in one job, see `launch.sh`.

## Outputs

Results are written to `--outdir/`:
- `<sample_id>_qcPlot.png`
- `<sample_id>_upsetPlot.png`
- `<sample_id>_metrics.txt`
- `<sample_id>_joint_qc.log` (multiome) or `<sample_id>_rna_qc.log` (scRNA/snRNA)

For `scRNA`/`snRNA` samples, `qcPlot.png` has fewer panels (no ATAC panels), and
`metrics.txt` has no `atac_*` columns. `scRNA` metrics additionally include
`rna_helm_metric`, `rna_log_helm_metric`, and `filter_helm` in place of
`filter_rna_exon_to_full_gene_body_ratio`.

## Container

The container is configured in `nextflow.config`:

```
process.container = 'library://alicewang24/python/jointqc_20260113:latest'
```

If you want to use a local `.sif`, replace the value with the absolute path and
ensure the directory is bound in `containerOptions`.
