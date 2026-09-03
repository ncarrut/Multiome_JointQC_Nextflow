# Multiome / scRNA / snRNA / ATAC Joint QC (Nextflow, cellranger branch)

This workflow submits 10X Cell Ranger (`cellranger-arc count`) data to a custom
QC analysis. Output includes background-corrected RNA-seq counts and a list of
pass-QC barcodes plus QC plots and visualizations.

Each sample's `location` always points at a `cellranger-arc count` outs/
directory (containing `raw_feature_bc_matrix.h5`, `gex_possorted_bam.bam`,
`atac_possorted_bam.bam`, `atac_fragments.tsv.gz`, etc). An optional `assay`
column in the samplesheet selects which QC module(s) that sample is routed
through:

- `multiome` full RNA+ATAC DAG, joint QC via
  `joint_qc.py`.
- `snRNA` — RNA half of the DAG only, QC via `rna_module_qc.py`
  (`--assay_res nucleus`), using the exon/full-gene-body intron-retention
  ratio filter.
- `scRNA` — RNA half of the DAG only, QC via `rna_module_qc.py`
  (`--assay_res cell`), using **HELM**
  (`rna_fraction_mitochondrial * (1 - rna_exon_to_full_gene_body_ratio)`) in
  place of the exon-ratio filter to capture the concordance between high exon-
  ratio and low mitochondrial reads. 
- `ATAC` — ATAC half of the DAG only, QC via `atac_module_qc.py`.

A single samplesheet can mix all four assay types — each sample is routed to
the matching set of processes automatically. `snRNA`/`scRNA` samples skip the
ATAC-side processes (`BIGWIG`, `PLOT_SIGNAL_AT_TSS`, `ATAQV_*`) and `JOINT_QC`;
`ATAC` samples skip the RNA-side processes (`SPLITTER`, `INTRONCOUNTER`, `QC`,
`CELLBENDER`, `EMPTYDROPS`) and the supplementary bigwig/bulk-ataqv processes,
running only `ATAQV_SINGLE_NUCLEUS` and `add_qc_metrics`.

## Requirements

- Nextflow (DSL2)
- Singularity/Apptainer module available on compute nodes

## Inputs

Required parameters:
- `--samplesheet` : tab-delimited text file with columns `sample` and
  `location` (path to the `cellranger-arc count` outs/ directory), 
   `assay` (`multiome`/`snRNA`/`scRNA`/`ATAC`) plus optionally
  `cluster_res`, and `df_pk`.
- `--results` : output directory for results (default: `results`).
- `--filter_MT_ATAC` : whether to filter ATAC nuclei on a %chrMT threshold;
  multiome and ATAC only (default: `true`).
- `--cellbender_fpr` : FPR threshold used for cellbender analysis
  (`<sample>.cellbender_FPR_<fpr>.h5`); multiome and RNA-only assays
  (default: `"0.05"`).

```
sample	assay	location	genome  cluster_res	df_pk
S1	multiome	/path/to/cellranger-arc/outs  hg38
S2	snRNA	/path/to/cellranger-arc/outs  hg38
S3	scRNA	/path/to/cellranger-arc/outs  hg38
S4	ATAC	/path/to/cellranger-arc/outs  hg38
```

## Run

```bash
nextflow run -resume path_to_main.nf \
  --samplesheet path_to_samplesheet.txt \
  --results results
```

## Output

`atac_qc`            -: ATAQV visualizations at bulk and single cell levels; `atac_qc/module_qc` has `atac_module_qc.py` outputs for `ATAC`-assay samples 

`atacv`              -: figures and data related to ATAQV 

`ataqv/single-nucleus` -: per-barcode ATAC QC metrics (`<sample>.txt`), input to `atac_module_qc.py`/`joint_qc.py` 

`bigwig`             -: ATAC bigwig files for TSS for selected genes (multiome only) 

`cellbender`         -: background filtered RNA counts: `<sample>.cellbender_FPR_0.05.h5`.  Use this for downstream analysis 

`counter`            -: intron/exon read counts (`<sample>_counts.txt`) used to derive the RNA exon/full-gene-body ratio 

`emptyDrops`         -: artifacts from Empty Drops 

`interactive-barcode...`     -: RNASeq barcode-rank plots (multiome only) 

`joint_qc`           -: plots and tables from the joint QC script for `multiome` samples. Pass-QC barcodes can be extracted from `<sample>_metrics.txt`. Start with `<sample>_qcPlot.png` to assess data quality; plots described below

`rna_qc`             -: `rna_module_qc.py` outputs for `snRNA`/`scRNA` samples 

`qc`                 -: RNASeq QC artifacts 

`splitter`           -: Separate RNA and ATAC data matrices. Note both modalities use the 'RNA barcodes' 

### Plot descriptions for `joint_qc/<sample>_qcPlot.png` 
**Top Row** 
* knee plot 
* mitochondrial fraction vs UMI 
* CellBender data, fraction of counts removed as ambient RNA vs UMI 

**Middle Row**
* **CellBender data**, number of droplets passing '% ambient removed' filter 
* **Cellbender data**, histogram of droplet cell probability for droplets passing EmptyDrops and mitochondrial fraction thresholds 
* Exon/full-gene-body count ratio vs UMI used to filter out droplets with too many exon vs intron counts 

**Bottom Row** 
* **ATAC data**, ATAC reads vs RNA UMIs 
* **ATAC data**, Transcription Start Site enrichment vs ATAC reads 
* **ATAC data**, ATAC mitochondrial read fraction vs ATAC reads.  This won't be used if `--filter_MT_ATAC` is set to false 


For `multiome`, `snRNA`, and `scRNA` samples, each QC module writes
`<sample>.qcPlot.png`, `<sample>.upsetPlot.png`, `<sample>.outmetrics.csv`
(or `<sample>_metrics.txt` for `joint_qc.py`), and `<sample>.log`.
`scRNA`/`snRNA` outmetrics have no `atac_*` columns; `ATAC` outmetrics have no
`rna_*` columns. `scRNA` metrics additionally include `rna_helm_metric`,
`rna_log_helm_metric`, and `filter_helm` in place of
`filter_rna_exon_to_full_gene_body_ratio`.
