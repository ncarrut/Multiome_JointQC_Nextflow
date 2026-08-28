# Multiome / scRNA / snRNA / ATAC Joint QC (Nextflow)

This workflow runs QC for four assay types, selected per-sample via an `assay`
column in the tab-delimited samplesheet:

- `multiome` (default if `assay` is omitted) — joint RNA/ATAC QC via `joint_qc.py`.
  Requires both `RNA_results_dir` and `ATAC_results_dir`.
- `snRNA` — RNA-only QC via `rna_module_qc.py` (`--assay_res nucleus`), using the
  standard exon/full-gene-body intron-retention ratio filter (same logic as
  multiome). Requires `RNA_results_dir`.
- `scRNA` — RNA-only QC via `rna_module_qc.py` (`--assay_res cell`), using **HELM**
  (`rna_fraction_mitochondrial * (1 - rna_exon_to_full_gene_body_ratio)`) in place
  of the exon-ratio filter, since whole-cell preps don't carry the same nuclear
  pre-mRNA contamination signature that the exon-ratio filter targets. Requires
  `RNA_results_dir`.
- `ATAC` — ATAC-only QC via `atac_module_qc.py`. Requires `ATAC_results_dir`.

A single library-config.tsv can mix all four assay types in one run — each
sample is routed to the matching process automatically. However, when the same
donor was profiled more than once (e.g. as both scRNA and snRNA, or across
different chemistry versions), each of those runs needs its own samplesheet and
`--outdir`, since the sample ID repeats and would otherwise overwrite outputs
from the other run. See `library-config_single_cell_v2.tsv`,
`library-config_single_cell_v3.tsv`, `library-config_single_nuclei_v2.tsv`,
`library-config_single_nuclei_v3.tsv`, and `library-config_single_nuclei_mixed.tsv`
for this project's actual per-type/chemistry samplesheets, and `launch.sh` for how
they're each run into a separate `results/<config>/` subdirectory in one job.

## Requirements

- Nextflow (DSL2)
- Singularity/Apptainer module available on compute nodes
- Container image (default: `docker://ncarrut/singlecell_qc:second`)

## Inputs

Required parameters:
- `--samplesheet` : tab-delimited file listing sample RNA (and, for multiome, ATAC) data locations for each sample (see `library-config_single_cell_v2.tsv` for an example)
- `--outdir` : output directory for results (default: `results`)
- `--filter_MT_ATAC` : whether to filter ATAC nuclei on a %chrMT threshold; multiome and ATAC only (default: `false`)
- `--cellbender_fpr` : FPR threshold used to locate the CellBender output file (`<sample>.cellbender_FPR_<fpr>.h5`); multiome and RNA-only assays (default: `0.05`)

The samplesheet is tab-delimited with a header row and these columns:
- `sample_id` : the library/sample ID.
- `assay` : one of `multiome`, `snRNA`, `scRNA`, `ATAC`. Optional, leave blank to default to `multiome`.
- `RNA_results_dir` : path to the upstream snRNAseq-NextFlow results directory. Required unless `assay` is `ATAC`.
- `ATAC_results_dir` : path to the upstream snATACseq-NextFlow results directory. Required when `assay` is `multiome` or `ATAC`.

```
sample_id	assay	RNA_results_dir	ATAC_results_dir
Ctrl-1_Pre	multiome	/path/to/snRNAseq-NextFlow/results/	/path/to/snATACseq-NextFlow/results/
SampleB	snRNA	/path/to/snRNAseq-NextFlow/results/
SampleC	scRNA	/path/to/snRNAseq-NextFlow/results/
SampleD	ATAC		/path/to/snATACseq-NextFlow/results/
```

## Run

```bash
nextflow run main.nf \
  --samplesheet library-config_single_cell_v2.tsv \
  --outdir results/single_cell_v2
```

To run every config for this project in one job, see `launch.sh`.

## Outputs

Results are written to `--outdir/`:
- `<sample_id>.qcPlot.png`
- `<sample_id>.upsetPlot.png`
- `<sample_id>.outmetrics.csv`
- `<sample_id>.log`

For `scRNA`/`snRNA` samples, `qcPlot.png` has fewer panels (no ATAC panels), and
`outmetrics.csv` has no `atac_*` columns. `scRNA` metrics additionally include
`rna_helm_metric`, `rna_log_helm_metric`, and `filter_helm` in place of
`filter_rna_exon_to_full_gene_body_ratio`. `ATAC` samples only have `atac_*`
columns and no RNA-side panels/columns.

## Container

The container is configured in `nextflow.config`:

```
process.container = 'library://alicewang24/python/jointqc_20260113:latest'
```

If you want to use a local `.sif`, replace the value with the absolute path and
ensure the directory is bound in `containerOptions`.
