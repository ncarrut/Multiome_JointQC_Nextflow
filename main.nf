#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.sample_id = null
params.rna_results_dir = null
params.atac_results_dir = null
params.outdir = "results"

process JOINT_QC {
    tag "${sample_id}"
    memory '32 GB'
    publishDir "${params.outdir}", mode: 'copy'

    input:
        tuple val(sample_id), val(rna_results_dir), val(atac_results_dir)

    output:
        path "${sample_id}_qcPlot.png", emit: qc_plot
        path "${sample_id}_upsetPlot.png", emit: upset_plot
        path "${sample_id}_metrics.txt", emit: metrics

    script:
    """
    python ${baseDir}/bin/joint_qc.py \\
        --sample ${sample_id} \\
        --RNA_results_dir ${rna_results_dir} \\
        --ATAC_results_dir ${atac_results_dir} \\
        --RNA_BARCODE_WHITELIST ${baseDir}/737K-arc-v1-rna.txt \\
        --ATAC_BARCODE_WHITELIST ${baseDir}/737K-arc-v1-atac.txt \\
        --qcPlot ${sample_id}_qcPlot.png \\
        --upsetPlot ${sample_id}_upsetPlot.png \\
        --outmetrics ${sample_id}_metrics.txt
    """
}

workflow {
    libraries = params.libraries.keySet()

    qc_pipeline_in = []

    for (library in libraries) {
        rna_results_dir = params.libraries[library]["RNA_results_dir"]
        atac_results_dir = params.libraries[library]["ATAC_results_dir"]
        qc_pipeline_in << [library, rna_results_dir, atac_results_dir]
    }

    qc_pipeline_ch = Channel.from(qc_pipeline_in)
    JOINT_QC(qc_pipeline_ch)
}
