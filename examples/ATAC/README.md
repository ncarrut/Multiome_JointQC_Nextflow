# Example: single-nucleus ATAC-seq

QC of single-nucleus chromatin accessibility only. ATAC-only samples
skip the RNA steps, CellBender, the bigwig steps and bulk ataqv.

## Input

Point `location` at a Cell Ranger output directory that contains:

| File | Used for |
|---|---|
| `atac_possorted_bam.bam` | ataqv per-nucleus ATAC metrics |

> **Note:** the pipeline expects the `cellranger-arc` BAM name
> (`atac_possorted_bam.bam`). Plain `cellranger-atac count` writes
> `possorted_bam.bam`. Symlink it to `atac_possorted_bam.bam` in the `outs/`
> directory, or in a directory of links that you use as `location`.

The `genome` column needs a `tss` entry in your config (see
[conf/references.example.config](../../conf/references.example.config)).

## Run

```bash
nextflow run /path/to/Multiome_JointQC_Nextflow/main.nf \
  --samplesheet examples/ATAC/samplesheet.txt \
  -c my_references.config -profile slurm

# Great Lakes
sbatch /path/to/Multiome_JointQC_Nextflow/conf/greatlakes/launch.sh samplesheet.txt
```

## Main outputs

| Directory | Contents |
|---|---|
| `atac_qc/module_qc/` | `<sample>.qcPlot.png`, `<sample>.upsetPlot.png`, `<sample>.outmetrics.csv` (per-barcode metrics and filter results), `<sample>.log` |
| `atac_qc/single_nucleus/` | Raw ataqv output |
| `ataqv/single-nucleus/` | Per-barcode ATAC metrics table (`<sample>.txt`) |
| `qc_report.html` | Summary report covering all samples |

To get the pass-QC barcodes, keep the rows of `atac_qc/module_qc/<sample>.outmetrics.csv`
where `pass_all_filters` is `True`.

## QC plot: `atac_qc/module_qc/<sample>.qcPlot.png`

<!-- Add an example image as examples/ATAC/qcPlot.png -->
![Example ATAC QC plot](qcPlot.png)

A 2×2 grid. Each point is one barcode, and dashed lines show the thresholds.
The title gives the number of pass-QC nuclei.

| Panel | Shows | Threshold |
|---|---|---|
| Top left | **TSS enrichment vs. HQAA** (high-quality autosomal alignments) | Red: min HQAA (multi-Otsu on log HQAA). Blue: min TSS enrichment = 2. |
| Top right | **ATAC barcode-rank plot**: HQAA per barcode, ranked | Red: min HQAA |
| Bottom left | **Mitochondrial % vs. HQAA** | Blue: max %chrMT. Drawn and applied only when `--filter_MT_ATAC true`. |
| Bottom right | **Max fraction of reads from a single autosome vs. HQAA**. A barcode whose reads come mostly from one chromosome is usually an artifact. | Blue: max fraction |

## Upset plot: `atac_qc/module_qc/<sample>.upsetPlot.png`

Counts how many barcodes pass each combination of filters. The bar with every
dot filled is the pass-QC set.

| Filter column | A barcode passes when… |
|---|---|
| `filter_atac_min_hqaa` | It has enough high-quality ATAC reads |
| `filter_atac_min_tss_enrichment` | TSS enrichment ≥ 2 |
| `filter_max_fraction_reads_from_single_autosome` | Reads aren't concentrated on one autosome |
| `filter_atac_max_mito` | Mitochondrial % is low (only with `--filter_MT_ATAC true`) |
