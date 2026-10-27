# Example: single-nucleus RNA-seq (snRNA)

QC of 10x 3' gene-expression libraries prepared from nuclei. A nucleus passes
when its RNA profile looks like an intact, uncontaminated nucleus.

## Input

Run the upstream snRNAseq-NextFlow pipeline first, then point `RNA_results_dir`
at its results directory. Leave `ATAC_results_dir` blank. The pipeline reads
these files (`<sample>` is the samplesheet's `sample_id`, which must match the
upstream file names):

| File | Used for |
|---|---|
| `cellbender/<sample>.cellbender_FPR_<fpr>.h5` | CellBender ambient fraction and cell probability (`<fpr>` = `--cellbender_fpr`) |
| `qc/<sample>.qc.txt` | Per-barcode RNA UMI and mitochondrial metrics |
| `starsolo/<sample>/<sample>.Solo.out/Gene/raw/` | Exon counts |
| `starsolo/<sample>/<sample>.Solo.out/GeneFull_ExonOverIntron/raw/` | Full-gene-body counts (exon ratio) |
| `emptyDrops/<sample>.knee.txt`, `<sample>.pass.txt` | Knee points and emptyDrops pass list |

> **Note:** end `RNA_results_dir` with a `/`. Some paths are built by plain
> string concatenation.

## Run

```bash
nextflow run /path/to/Multiome_JointQC_Nextflow/main.nf \
  --samplesheet examples/single_nuclei_RNASeq/samplesheet.txt \
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
where `pass_all_filters` is `True`. For downstream analysis, use the upstream
CellBender counts (`<RNA_results_dir>/cellbender/`).

## QC plot: `<sample>.qcPlot.png`

![Example snRNA QC plot](single_nuclei_example_qcPlot.png)

A 2×3 grid; it's the RNA part of the multiome plot. The title gives the number
of pass-QC nuclei.

| Panel | Shows | Threshold |
|---|---|---|
| Top left | **Barcode-rank (knee) plot**: RNA UMIs per barcode, ranked | Knee, inflection, end-of-cliff and plateau lines from emptyDrops. The title shows the number of knees found. |
| Top middle | **RNA mitochondrial fraction vs. UMIs** | Blue: max mito fraction. Red: min UMIs. |
| Top right | **CellBender fraction of counts removed vs. UMIs** | Blue: max fraction removed |
| Bottom left | **Histogram of % ambient removed** (barcodes between 5% and 50%) | Blue: multi-Otsu cutoff |
| Bottom middle | **CellBender cell probability** for barcodes that pass emptyDrops and the mito filter | Barcodes need probability ≥ 0.99 |
| Bottom right | **Exon / full-gene-body ratio vs. UMIs** | Red: max ratio. Nuclear RNA is intron-rich, so barcodes with mostly exonic reads are likely ambient or cytoplasmic contamination. |

## Upset plot: `<sample>.upsetPlot.png`

Counts how many barcodes pass each combination of filters. The bar with every
dot filled is the pass-QC set.

| Filter column | A barcode passes when… |
|---|---|
| `filter_rna_emptyDrops` | The upstream emptyDrops run distinguishes it from background |
| `filter_rna_min_umi` | It has enough RNA UMIs |
| `filter_pct_cellbender_removed` | Ambient RNA removed by CellBender is below the cutoff |
| `filter_cellbender_cell_probability` | CellBender cell probability ≥ 0.99 |
| `filter_rna_max_mito` | RNA mitochondrial fraction is low |
| `filter_rna_exon_to_full_gene_body_ratio` | Exon/full-gene-body ratio is below the cutoff |
