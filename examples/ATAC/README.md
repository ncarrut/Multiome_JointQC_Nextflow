# Example: single-nucleus ATAC-seq

QC of single-nucleus chromatin accessibility only.

## Input

Run the upstream snATACseq-NextFlow pipeline first, then point `ATAC_results_dir`
at its results directory. Leave `RNA_results_dir` blank. The pipeline reads:

| File | Used for |
|---|---|
| `ataqv/single-nucleus/<sample>.txt` | ataqv per-nucleus ATAC metrics |

`<sample>` is the samplesheet's `sample_id` and must match the upstream file name.

## Run

```bash
nextflow run /path/to/Multiome_JointQC_Nextflow/main.nf \
  --samplesheet examples/ATAC/samplesheet.txt \
  --outdir results -profile slurm

# Great Lakes
sbatch /path/to/Multiome_JointQC_Nextflow/conf/greatlakes/launch.sh samplesheet.txt [outdir]
```

## Main outputs

All files are written directly to `--outdir`:

| File | Contents |
|---|---|
| `<sample>.qcPlot.png` | QC plot (below) |
| `<sample>.upsetPlot.png` | Upset plot (below) |
| `<sample>.outmetrics.csv` | Per-barcode metrics and filter results |
| `<sample>.log` | Run log |

To get the pass-QC barcodes, keep the rows of `<sample>.outmetrics.csv`
where `pass_all_filters` is `True`.

## QC plot: `<sample>.qcPlot.png`

<!-- Add an example image as examples/ATAC/qcPlot.png and link it here -->

A 2×2 grid. Each point is one barcode, and dashed lines show the thresholds.
The title gives the number of pass-QC nuclei.

| Panel | Shows | Threshold |
|---|---|---|
| Top left | **TSS enrichment vs. HQAA** (high-quality autosomal alignments) | Red: min HQAA (multi-Otsu on log HQAA). Blue: min TSS enrichment = 2. |
| Top right | **ATAC barcode-rank plot**: HQAA per barcode, ranked | Red: min HQAA |
| Bottom left | **Mitochondrial % vs. HQAA** | Blue: max %chrMT. Drawn and applied only when `--filter_MT_ATAC true`. |
| Bottom right | **Max fraction of reads from a single autosome vs. HQAA**. A barcode whose reads come mostly from one chromosome is usually an artifact. | Blue: max fraction |

## Upset plot: `<sample>.upsetPlot.png`

Counts how many barcodes pass each combination of filters. The bar with every
dot filled is the pass-QC set.

| Filter column | A barcode passes when… |
|---|---|
| `filter_atac_min_hqaa` | It has enough high-quality ATAC reads |
| `filter_atac_min_tss_enrichment` | TSS enrichment ≥ 2 |
| `filter_max_fraction_reads_from_single_autosome` | Reads aren't concentrated on one autosome |
| `filter_atac_max_mito` | Mitochondrial % is low (only with `--filter_MT_ATAC true`) |
