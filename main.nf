#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.sample_id = null
params.rna_results_dir = null
params.atac_results_dir = null
params.outdir = "results"

process JOINT_QC {
    tag "${sample_id}"
    memory '32 GB'
    container 'docker://porchard/general:20220406125608'
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
        tuple val(sample_id), val(rna_results_dir), val(atac_results_dir)

    output:
        path "qcPlot.png", emit: qc_plot
        path "upsetPlot.png", emit: upset_plot
        path "metrics.txt", emit: metrics

    script:
    """
    python ${baseDir}/bin/joint_qc.py \\
        --sample ${sample_id} \\
        --RNA_results_dir ${rna_results_dir} \\
        --ATAC_results_dir ${atac_results_dir} \\
        --RNA_BARCODE_WHITELIST ${baseDir}/737K-arc-v1-rna.txt \\
        --ATAC_BARCODE_WHITELIST ${baseDir}/737K-arc-v1-atac.txt \\
        --qcPlot qcPlot.png \\
        --upsetPlot upsetPlot.png \\
        --outmetrics metrics.txt
    """
}

workflow {
    if (!params.sample_id || !params.rna_results_dir || !params.atac_results_dir) {
        error "Missing params. Example: nextflow run main.nf --sample_id S1 --rna_results_dir /path/RNA --atac_results_dir /path/ATAC"
    }

    JOINT_QC(tuple(params.sample_id, params.rna_results_dir, params.atac_results_dir))
}
