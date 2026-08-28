#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process JOINT_QC {
    tag "${sample_id}"
    memory '32 GB'
    publishDir "${params.outdir}", mode: 'copy'

    input:
        tuple val(sample_id), val(rna_results_dir), val(atac_results_dir)

    output:
        path "${sample_id}.qcPlot.png", emit: qc_plot
        path "${sample_id}.upsetPlot.png", emit: upset_plot
        path "${sample_id}.outmetrics.csv", emit: metrics
        path "${sample_id}.log", emit: log

    script:
    """
    python ${baseDir}/bin/joint_qc.py \\
        --sample ${sample_id} \\
        --RNA_results_dir ${rna_results_dir} \\
        --cellbender_fpr ${params.cellbender_fpr} \\
        --ATAC_results_dir ${atac_results_dir} \\
        --RNA_BARCODE_WHITELIST ${baseDir}/737K-arc-v1-rna.txt \\
        --ATAC_BARCODE_WHITELIST ${baseDir}/737K-arc-v1-atac.txt \\
        --filter_MT_ATAC ${params.filter_MT_ATAC} \\
        --output .
    """
}

process RNA_QC {
    tag "${sample_id}"
    memory '16 GB'
    publishDir "${params.outdir}", mode: 'copy'

    input:
        tuple val(sample_id), val(rna_results_dir), val(assay_res)

    output:
        path "${sample_id}.qcPlot.png", emit: qc_plot
        path "${sample_id}.upsetPlot.png", emit: upset_plot
        path "${sample_id}.outmetrics.csv", emit: metrics
        path "${sample_id}.log", emit: log

    script:
    """
    python ${baseDir}/bin/rna_module_qc.py \\
        --sample ${sample_id} \\
        --assay_res ${assay_res} \\
        --RNA_results_dir ${rna_results_dir} \\
        --cellbender_fpr ${params.cellbender_fpr} \\
        --output .
    """
}

process ATAC_QC {
    tag "${sample_id}"
    memory '16 GB'
    publishDir "${params.outdir}", mode: 'copy'

    input:
        tuple val(sample_id), val(atac_results_dir)

    output:
        path "${sample_id}.qcPlot.png", emit: qc_plot
        path "${sample_id}.upsetPlot.png", emit: upset_plot
        path "${sample_id}.outmetrics.csv", emit: metrics
        path "${sample_id}.log", emit: log

    script:
    """
    python ${baseDir}/bin/atac_module_qc.py \\
        --sample ${sample_id} \\
        --ATAC_results_dir ${atac_results_dir} \\
        --filter_MT_ATAC ${params.filter_MT_ATAC} \\
        --output .
    """
}

workflow {
    def valid_assays = ["multiome", "scRNA", "snRNA", "ATAC"]
    def assay_res = ["scRNA": "cell", "snRNA": "nucleus"]

    multiome_in = []
    rna_in = []
    atac_in = []

    def rows = file(params.samplesheet).splitCsv(header: true, sep: '\t')

    for (row in rows) {
        library = row.sample_id
        assay = row.assay ? row.assay : "multiome"

        if (!valid_assays.contains(assay)) {
            error "Library '${library}' has unrecognized assay '${assay}'. Expected one of: ${valid_assays.join(', ')}."
        }

        if (assay == "ATAC") {
            atac_results_dir = row.ATAC_results_dir
            if (!atac_results_dir) {
                error "Library '${library}' has assay 'ATAC' but is missing ATAC_results_dir."
            }
            atac_in << [library, atac_results_dir]
            continue
        }

        rna_results_dir = row.RNA_results_dir
        if (!rna_results_dir) {
            error "Library '${library}' is missing RNA_results_dir."
        }

        if (assay == "multiome") {
            atac_results_dir = row.ATAC_results_dir
            if (!atac_results_dir) {
                error "Library '${library}' has assay 'multiome' but is missing ATAC_results_dir."
            }
            multiome_in << [library, rna_results_dir, atac_results_dir]
        } else {
            rna_in << [library, rna_results_dir, assay_res[assay]]
        }
    }

    if (multiome_in) {
        multiome_ch = Channel.from(multiome_in)
        JOINT_QC(multiome_ch)
    }

    if (rna_in) {
        rna_ch = Channel.from(rna_in)
        RNA_QC(rna_ch)
    }

    if (atac_in) {
        atac_ch = Channel.from(atac_in)
        ATAC_QC(atac_ch)
    }
}
