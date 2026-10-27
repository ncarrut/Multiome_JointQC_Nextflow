# Multiome / scRNA / snRNA / ATAC Joint QC (Nextflow)

A Nextflow pipeline that runs barcode-level QC on the results of the upstream
snRNAseq-NextFlow and snATACseq-NextFlow pipelines. It works with
**Multiome (RNA + ATAC)**, **single-nucleus RNA**, **single-cell RNA** and
**single-nucleus ATAC** data. For each sample it produces:

- a per-barcode metrics table with a `pass_all_filters` column
- diagnostic QC and upset plots

Thresholds (UMI, mitochondrial fraction, ambient fraction, exon ratio, ATAC
depth, …) are set automatically from each sample's own distributions. None
need to be tuned by hand.

## Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2, ≥ 20.10)
- Singularity or Apptainer
- Upstream results: snRNAseq-NextFlow (with CellBender and emptyDrops) for RNA,
  snATACseq-NextFlow (with ataqv) for ATAC

## Quick start

1. Pick the example that matches your data and copy its samplesheet:

   | Data type | Example |
   |---|---|
   | 10x Multiome | [examples/single_nuclei_multiome](examples/single_nuclei_multiome/README.md) |
   | Single-nucleus RNA | [examples/single_nuclei_RNASeq](examples/single_nuclei_RNASeq/README.md) |
   | Single-cell RNA | [examples/single_cell_RNASeq](examples/single_cell_RNASeq/README.md) |
   | Single-nucleus ATAC | [examples/ATAC](examples/ATAC/README.md) |
   | Several types in one run | [examples/mixed/samplesheet.txt](examples/mixed/samplesheet.txt) |

   Each example README lists the upstream files needed and explains every
   panel of the QC plots.

2. Run:

   ```bash
   nextflow run /path/to/Multiome_JointQC_Nextflow/main.nf \
     --samplesheet samplesheet.txt \
     --outdir results \
     -profile slurm
   ```

## Samplesheet

Tab-delimited with a header row, one row per sample.

| Column | Required | Description |
|---|---|---|
| `sample_id` | yes | Sample name. Must match the upstream file names; also used in output file names. |
| `assay` | no (default `multiome`) | `multiome`, `snRNA`, `scRNA` or `ATAC` |
| `RNA_results_dir` | unless `ATAC` | snRNAseq-NextFlow results directory (end it with `/`) |
| `ATAC_results_dir` | `multiome`, `ATAC` | snATACseq-NextFlow results directory |

```
sample_id	assay	RNA_results_dir	ATAC_results_dir
M1	multiome	/path/to/snRNAseq-NextFlow/results/	/path/to/snATACseq-NextFlow/results/
N1	snRNA	/path/to/snRNAseq-NextFlow/results/
C1	scRNA	/path/to/snRNAseq-NextFlow/results/
A1	ATAC		/path/to/snATACseq-NextFlow/results/
```

One samplesheet can mix all four assay types. However, if the same donor was
profiled more than once (e.g. as both scRNA and snRNA, or with different
chemistry versions), give each run its own samplesheet and `--outdir`.
Otherwise the repeated sample ID overwrites the other run's outputs.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--samplesheet` | — | Path to the samplesheet (required) |
| `--outdir` | `results` | Output directory |
| `--filter_MT_ATAC` | `false` | Filter nuclei on ATAC mitochondrial % (multiome and ATAC) |
| `--cellbender_fpr` | `0.05` | FPR of the upstream CellBender run, used to find `cellbender/<sample>.cellbender_FPR_<fpr>.h5` (multiome and RNA) |

## Profiles and site configuration

| Profile | Use |
|---|---|
| `standard` (default) | Run every process on the local machine |
| `slurm` | Submit each process as a SLURM job. Add account/partition options in your own config passed with `-c`. |
| `greatlakes` | University of Michigan Great Lakes cluster (account, bind mounts) |

Site-specific settings go in a separate config file. See
[conf/greatlakes.config](conf/greatlakes.config) for a complete example.

### Great Lakes users

Submit the driver job from the directory where the output directory and `work/` should be
written:

```bash
sbatch /path/to/Multiome_JointQC_Nextflow/conf/greatlakes/launch.sh samplesheet.txt [outdir]
```

## How samples are processed

Each sample is routed by its `assay`:

| Assay | QC script | Notes |
|---|---|---|
| `multiome` | `joint_qc.py` | Joint RNA + ATAC QC; RNA and ATAC barcodes matched with the bundled 737K-arc-v1 whitelists |
| `snRNA` | `rna_module_qc.py --assay_res nucleus` | Filters on exon / full-gene-body ratio |
| `scRNA` | `rna_module_qc.py --assay_res cell` | Filters on HELM (`rna_fraction_mitochondrial × (1 − rna_exon_to_full_gene_body_ratio)`) instead of exon ratio, since whole cells are naturally exon-rich |
| `ATAC` | `atac_module_qc.py` | ataqv per-nucleus metrics only |

## Outputs

All files are written directly to `--outdir`:

| File | Contents |
|---|---|
| `<sample_id>.qcPlot.png` | Diagnostic QC plot |
| `<sample_id>.upsetPlot.png` | Barcodes passing each combination of filters |
| `<sample_id>.outmetrics.csv` | Per-barcode metrics and filter results |
| `<sample_id>.log` | Run log |

The pass-QC barcodes are the rows where `pass_all_filters` is `True`.
`scRNA`/`snRNA` plots have no ATAC panels and their metrics have no `atac_*`
columns. `scRNA` metrics have `rna_helm_metric`, `rna_log_helm_metric` and
`filter_helm` in place of `filter_rna_exon_to_full_gene_body_ratio`. `ATAC`
samples have only `atac_*` columns. For an explanation of the plots, see the
example README for your data type.

## Container

The container is set in `nextflow.config`:

```
process.container = 'docker://ncarrut/singlecell_qc:second'
```

To use a local `.sif`, replace the value with its absolute path and make sure
its directory is bound in `containerOptions`.
