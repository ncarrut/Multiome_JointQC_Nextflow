#!/usr/bin/env nextflow

nextflow.enable.dsl=2

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
        path "${sample_id}_joint_qc.log", emit: log

    script:
    """
    python ${baseDir}/bin/joint_qc.py \\
        --sample ${sample_id} \\
        --RNA_results_dir ${rna_results_dir} \\
        --ATAC_results_dir ${atac_results_dir} \\
        --RNA_BARCODE_WHITELIST ${baseDir}/737K-arc-v1-rna.txt \\
        --ATAC_BARCODE_WHITELIST ${baseDir}/737K-arc-v1-atac.txt \\
        --filter_MT_ATAC ${params.filter_MT_ATAC} \\
        --qcPlot ${sample_id}_qcPlot.png \\
        --upsetPlot ${sample_id}_upsetPlot.png \\
        --outmetrics ${sample_id}_metrics.txt \\
        --outlogs ${sample_id}_joint_qc.log
    """
}

process RNA_QC {
    tag "${sample_id}"
    memory '16 GB'
    publishDir "${params.outdir}", mode: 'copy'

    input:
        tuple val(sample_id), val(rna_results_dir), val(assay)

    output:
        path "${sample_id}_qcPlot.png", emit: qc_plot
        path "${sample_id}_upsetPlot.png", emit: upset_plot
        path "${sample_id}_metrics.txt", emit: metrics
        path "${sample_id}_rna_qc.log", emit: log

    script:
    """
    python ${baseDir}/bin/rna_qc.py \\
        --sample ${sample_id} \\
        --RNA_results_dir ${rna_results_dir} \\
        --assay ${assay} \\
        --qcPlot ${sample_id}_qcPlot.png \\
        --upsetPlot ${sample_id}_upsetPlot.png \\
        --outmetrics ${sample_id}_metrics.txt \\
        --outlogs ${sample_id}_rna_qc.log
    """
}

workflow {
    libraries = params.libraries.keySet()

    def valid_assays = ["multiome", "scRNA", "snRNA"]

    multiome_in = []
    rna_in = []

    for (library in libraries) {
        entry = params.libraries[library]
        assay = entry.containsKey("assay") ? entry["assay"] : "multiome"

        if (!valid_assays.contains(assay)) {
            error "Library '${library}' has unrecognized assay '${assay}'. Expected one of: ${valid_assays.join(', ')}."
        }

        rna_results_dir = entry["RNA_results_dir"]
        if (!rna_results_dir) {
            error "Library '${library}' is missing RNA_results_dir."
        }

        if (assay == "multiome") {
            atac_results_dir = entry["ATAC_results_dir"]
            if (!atac_results_dir) {
                error "Library '${library}' has assay 'multiome' but is missing ATAC_results_dir."
            }
            multiome_in << [library, rna_results_dir, atac_results_dir]
        } else {
            rna_in << [library, rna_results_dir, assay]
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
}
