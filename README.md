# Multiome Joint QC (Nextflow)

This workflow submits 10X cellranger multiome data to a custom QC analysis. Output includes background corrected RNASeq counts and a list of passQC barcodes plus QC plots and visualizations.  

## Requirements

- Nextflow (DSL2)
- Singularity/Apptainer module available on compute nodes

## Inputs

Required parameters:
- `--samplesheet` : tab-delimited text file with at least the columns `sample` and `location` giving sample names and locations for the cellranger output.

## Run

```bash
nextflow run -resume path_to_main.nf \
  --samplesheet path_to_samplesheet.txt \
  --results results
```

## Output
`atac_qc`         -: ATAQV visualizations at bulk and single cell levels
`atacv`           -: figures and data related to ATAQV
`bigwig`          -: ATAC bigwig files for TSS for selected genes
`cellbender`      -: background filtered RNA counts:  `<sample>.cellbender_FPR_0.05.h5` 
`counter`         -: intermediate data for intron/exon quantification
`emptyDrops`      -: artifacts from Empty Drops
`interactive...`  -: RNASeq barcode-rank plots
`joint_qc`        -: plots and tables from the joint QC script.  Passqc barcodes can be extracted from `<sample>_metrics.txt`
`qc`              -: RNASeq QC artifacts 
`splitter`        -: Separate RNA and ATAC data matrices.  Note both modalities use the 'RNA barcodes'

