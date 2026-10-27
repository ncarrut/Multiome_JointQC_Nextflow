# Single-cell / single-nucleus Joint QC (Nextflow)

A Nextflow pipeline that runs barcode-level QC on 10x Genomics Cell Ranger
output. It works with **Multiome (RNA + ATAC)**, **single-nucleus RNA**,
**single-cell RNA** and **single-nucleus ATAC** data. For each sample it
produces:

- ambient-RNA-corrected counts (CellBender), for downstream analysis
- a per-barcode metrics table with a `pass_all_filters` column
- diagnostic QC and upset plots
- `qc_report.html`, a summary report covering all samples

Thresholds (UMI, mitochondrial fraction, ambient fraction, exon ratio, ATAC
depth, …) are set automatically from each sample's own distributions. None
need to be tuned by hand.

## Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2, ≥ 20.10)
- Singularity or Apptainer 
- A GPU node for the CellBender step

## Quick start

1. Pick the example that matches your data and copy its samplesheet:

   | Data type | Example |
   |---|---|
   | 10x Multiome (cellranger-arc) | [examples/single_nuclei_multiome](examples/single_nuclei_multiome/README.md) |
   | Single-nucleus RNA (cellranger) | [examples/single_nuclei_RNASeq](examples/single_nuclei_RNASeq/README.md) |
   | Single-cell RNA (cellranger) | [examples/single_cell_RNASeq](examples/single_cell_RNASeq/README.md) |
   | Single-nucleus ATAC | [examples/ATAC](examples/ATAC/README.md) |
   | Several types in one run | [examples/mixed/samplesheet.txt](examples/mixed/samplesheet.txt) |

   Each example README lists the input files needed and explains every panel
   of the QC plots.

2. For multiome or ATAC samples, copy
   [conf/references.example.config](conf/references.example.config) and fill in the
   TSS and chromosome-sizes files for your genome(s).

3. Run:

   ```bash
   nextflow run -resume /path/to/Multiome_JointQC_Nextflow/main.nf \
     --samplesheet samplesheet.txt \
     --results results \
     -c my_references.config \
     -profile slurm
   ```

## Samplesheet

Tab-delimited, one row per sample.

| Column | Required | Description |
|---|---|---|
| `sample` | yes | Sample name, used in output file names |
| `assay` | no (default `multiome`) | `multiome`, `snRNA`, `scRNA` or `ATAC` |
| `location` | yes | Cell Ranger `outs/` directory |
| `genome` | yes | Genome key, e.g. `hg38`, `mm10`. Supported: hg19, hg38, mm9, mm10, rn4–rn7. |
| `cluster_res` | no | Overrides `--cluster_res` for this sample |
| `df_pk` | no | Overrides `--df_pk` for this sample |

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--samplesheet` | — | Path to the samplesheet (required) |
| `--results` | `results` | Output directory |
| `--filter_MT_ATAC` | `true` | Filter nuclei on ATAC mitochondrial % (multiome and ATAC) |
| `--emptydrops_fdr` | `0.01` | FDR cutoff for calling a barcode non-empty with emptyDrops |
| `--cellbender_fpr` | `0.05` | Which CellBender FPR result to use downstream. Must be one of `params.cellbender_fprs` (`0.01`, `0.05`); CellBender writes all of them to `cellbender/`. |
| `params.tss` | `[:]` | Map of genome → TSS BED (multiome/ATAC). Set in a config file. |
| `params.chrom_sizes` | `[:]` | Map of genome → chromosome sizes (multiome). Set in a config file. |
| `params.plot_signal_at_genes` | GAPDH, MYH1, MYH2, MYH7 | Genes whose TSS ATAC signal is plotted (multiome) |

## Profiles and site configuration

| Profile | Use |
|---|---|
| `standard` (default) | Run every process on the local machine |
| `slurm` | Submit each process as a SLURM job. Add account/partition options in your own config passed with `-c`. |
| `greatlakes` | University of Michigan Great Lakes cluster (accounts, GPU partition, bind mounts, reference paths) |

Site-specific settings go in a separate config file. See
[conf/greatlakes.config](conf/greatlakes.config) for a complete example.

### Great Lakes users

The Great Lakes profile includes the reference files, so no `-c` is needed.
Submit the driver job from the directory where `results/` and `work/` should
be written:

```bash
sbatch /path/to/Multiome_JointQC_Nextflow/conf/greatlakes/launch.sh samplesheet.txt [results_dir]
```

## How samples are processed

Each sample is routed by its `assay`:

| Assay | Steps | QC script | Results |
|---|---|---|---|
| `multiome` | RNA + ATAC: CellBender, emptyDrops, intron counting, ataqv, bigwigs | `joint_qc.py` | `joint_qc/` |
| `snRNA` | RNA only; filters on exon / full-gene-body ratio | `rna_module_qc.py --assay_res nucleus` | `rna_qc/` |
| `scRNA` | RNA only; filters on HELM instead of exon ratio | `rna_module_qc.py --assay_res cell` | `rna_qc/` |
| `ATAC` | ataqv per-nucleus metrics only | `atac_module_qc.py` | `atac_qc/module_qc/` |

## Outputs

| Directory | Contents | Assays |
|---|---|---|
| `qc_report.html` | Summary report: filters applied, cells passing per sample, QC plots | all |
| `joint_qc/` | `<sample>_qcPlot.png`, `<sample>_upsetPlot.png`, `<sample>_metrics.txt` | multiome |
| `rna_qc/` | `<sample>.qcPlot.png`, `<sample>.upsetPlot.png`, `<sample>.outmetrics.csv` | snRNA, scRNA |
| `atac_qc/module_qc/` | `<sample>.qcPlot.png`, `<sample>.upsetPlot.png`, `<sample>.outmetrics.csv` | ATAC |
| `cellbender/` | Ambient-corrected counts `<sample>.cellbender_FPR_<fpr>.h5`. **Use for downstream analysis.** | multiome, RNA |
| `emptyDrops/` | emptyDrops knee points and pass list | multiome, RNA |
| `counter/` | Per-barcode intron/exon read counts | multiome, RNA |
| `qc/` | Per-barcode RNA read/UMI metrics | multiome, RNA |
| `splitter/` | RNA (GEX) and ATAC matrices split from the Cell Ranger h5 | multiome, RNA |
| `interactive-barcode-rank-plots/` | Interactive RNA barcode-rank plots | multiome |
| `atac_qc/`, `ataqv/single-nucleus/` | ataqv metrics: bulk HTML viewer, per-nucleus tables and plots | multiome, ATAC |
| `bigwig/` | ATAC coverage bigwigs and TSS signal plots | multiome |

In every metrics table, the pass-QC barcodes are the rows where
`pass_all_filters` is `True`. For an explanation of the plots, see the
example README for your data type.
