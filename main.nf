#!/usr/bin/env nextflow

nextflow.enable.dsl=2

ORGANISMS = ['hg19': 'human', 
    'hg38': 'human',
    'rn4': 'rat',
    'rn5': 'rat',
    'rn6': 'rat',
    'rn7': 'rat',
    'mm9': 'mouse',
    'mm10': 'mouse'
]

def get_organism (genome) {
    return ORGANISMS[genome]
}

def get_tss(genome) {
    return params.tss[genome] ?: params.tss['hg38']
}

def get_chrom_sizes(genome) {
    return params.chrom_sizes[genome] ?: params.chrom_sizes['hg38']
}

process SPLITTER{
    cpus 1
    memory '40 GB'
    publishDir "${params.results}/splitter/", mode: 'link', overwrite: true, saveAs: { filename ->
        if (filename.contains("_GEX")) "GEX/$filename"
        else if (filename.contains("_ATAC")) "ATAC/$filename"
        else if (filename.endsWith(".log")) "logs/$filename"
        else null
    }
    container 'docker://ncarrut/cellbender-scanpy:multivelo_updated'
    time '2h'
    errorStrategy 'ignore'
    
    input:
    tuple val(sample), val(location), val(cluster_res), val(df_pk)
    output:
    tuple val(sample), path("${sample}_GEX"), emit: GEX
    tuple val(sample), path("${sample}_ATAC"), emit: ATAC
    path("*.log")

    """
    splitter.py --sample ${sample} --h5 "${location}/raw_feature_bc_matrix.h5" 2>&1 | tee ${sample}.log
    """
}

process INTRONCOUNTER{
    cpus 1
    memory '32 GB'
    publishDir "${params.results}/counter/", mode: 'link', overwrite: true
    container 'docker://ncarrut/samtools:20260405'
    time '4h'
    errorStrategy 'ignore'

    input:
    tuple val(sample), val(location), val(cluster_res), val(df_pk)
    output:
    tuple val(sample), path("${sample}_counts.txt"), emit: counts
    path("*.log")

    """
    intron_counter.sh ${location}/gex_possorted_bam.bam ${sample}_counts.txt 2>&1 | tee ${sample}.log
    """
}

process QC{

    memory '25 GB'
    publishDir "${params.results}/qc"
    tag "${sample}"
    container 'library://porchard/default/general:20220107'
    cpus 1
    time '5h'

    input:
    tuple val(sample), val(bam_location), val(cluster_res), val(df_pk), val(mtx_location)
 
    output:
    tuple val(sample), path("${sample}.qc.txt")

    """
    qc-from-starsolo.py ${bam_location}/gex_possorted_bam.bam ${mtx_location}/matrix.mtx ${mtx_location}/barcodes.tsv > ${sample}.qc.txt
    """

}


process PLOTQC {

    memory '15 GB'
    publishDir "${params.results}/qc"
    tag "${sample}"
    container 'library://porchard/default/dropkick:20220225'
    cpus 1
    time '5h'

    input:
    tuple val(sample), path(metrics)

    output:
    tuple val(sample), path("${sample}metrics.png"), path("${sample}suggested-thresholds.tsv")


    """
    plot-qc-metrics.py --prefix ${sample} ${metrics}
    """

}


process INTERACTIVEBARCODERANKPLOT {

    memory '15 GB'
    publishDir "${params.results}/interactive-barcode-rank-plots"
    tag "${sample}"
    container "docker://porchard/plotly:20230705"
    cpus 1
    time '3h'

    input:
    tuple val(sample), val(location)

    output:
    path("${sample}.barcode-rank-plot.html")


    """
    interactive-barcode-rank-plot.py ${location}/matrix.mtx ${sample}.barcode-rank-plot.html
    """
}

process CELLBENDER{
    cpus 1
    memory '40 GB'
    publishDir "${params.results}/cellbender/", mode: 'link', overwrite: true
    container 'docker://porchard/cellbender:0.3.0'
    time '12h'
    errorStrategy 'ignore'

    input:
    tuple val(sample), val(location)
    output:
    path("${sample}.*")
    tuple val(sample), path("${sample}*.h5"), emit: h5_files
    tuple val(sample), path("${sample}.cellbender_FPR_0.05.h5"), emit: h5_fpr05
    tuple val(sample), path("${sample}.cellbender_FPR_0.05_metrics.csv"), emit: metrics
    path("*.log")

    """
    cellbender remove-background --cuda --epochs 150 --fpr 0.01 0.05 --input ${location} --output ./${sample}.cellbender.h5 \
        --projected-ambient-count-threshold 2 --posterior-batch-size 64 --exclude-feature-types Peaks
    cp .command.log ${sample}.log
    """
}

process EMPTYDROPS {
    publishDir "${params.results}/emptyDrops", mode: 'link', overwrite: true
    container 'docker://porchard/dropletutils:20241202'
    time '4h'
    errorStrategy 'ignore'
    memory '8 GB'

    input:
    tuple val(sample), val(location), path(cellbender_metrics)

    output:
    tuple val(sample), path("${sample}.knee.txt"), path("${sample}.pass.txt"), emit: knee_pass

    """
    emptyDrops_wCellBender.R --donor ${sample} --barcodeList ${location} --cbMetrics ${cellbender_metrics} --lowerForKnee 100 --outKnee ${sample}.knee.txt --outPass ${sample}.pass.txt
    """
}


process BIGWIG {

    time '24h'
    publishDir "${params.results}/bigwig"
    tag "${sample}"
    container 'library://porchard/default/general:20220107'
    memory { 20.GB * task.attempt }
    errorStrategy 'ignore'
    cpus 1
    label 'largemem'

    input:
    tuple val(sample), val(location), val(cluster_res), val(df_pk)

    output:
    tuple val(sample), path("${sample}.bw")

    """
    # Convert fragments to bedgraph for bigwig generation
    zcat ${location}/atac_fragments.tsv.gz | awk 'BEGIN{OFS="\t"} {print \$1,\$2,\$3,\$4}' | sort -k1,1 -k2,2n | bedtools merge -i - -c 4 -o count > ${sample}.bedgraph
    LC_COLLATE=C sort -k1,1 -k2n,2 ${sample}.bedgraph > sorted.bedgraph
    bedClip sorted.bedgraph ${get_chrom_sizes(params.genome)} clipped.bedgraph
    bedGraphToBigWig clipped.bedgraph ${get_chrom_sizes(params.genome)} ${sample}.bw
    rm sorted.bedgraph clipped.bedgraph ${sample}.bedgraph
    """

}


process PLOT_SIGNAL_AT_TSS {

    publishDir "${params.results}/bigwig/plot"
    errorStrategy 'retry'
    maxRetries 1
    memory { 10.GB * task.attempt }
    tag "${sample}"
    container 'library://porchard/default/general:20220107'
    cpus 1
    time '24h'

    input:
    tuple val(sample), path(bw)

    output:
    path("*.png") optional true

    """
    plot-signal-at-tss.py --genes ${params.plot_signal_at_genes.join(' ')} --tss-file ${get_tss(params.genome)} --bigwigs ${bw}
    """

}


process ATAQV_SINGLE_NUCLEUS {
    errorStrategy 'retry'
    maxRetries 1
    memory { 5.GB * task.attempt }
    time '12h'
    tag "${sample}"
    container 'docker://porchard/ataqv:1.5.0'
    cpus 2
    publishDir "${params.results}/atac_qc/single_nucleus", mode: 'link'

    input:
    tuple val(sample), val(location), val(cluster_res), val(df_pk)
    
    output:
    tuple val(sample), path("${sample}.ataqv.txt.gz"), emit: metrics
    path("${sample}.ataqv.out")

    script:
    """
    ataqv --name ${sample} --ignore-read-groups --nucleus-barcode-tag CB --metrics-file ${sample}.ataqv.txt.gz --tss-file ${get_tss(params.genome)} ${get_organism(params.genome)} ${location}/atac_possorted_bam.bam > ${sample}.ataqv.out
    """
}

process add_qc_metrics {

    publishDir "${params.results}/ataqv/single-nucleus"
    time '1h'
    tag "${library}"
    container 'library://porchard/default/general:20220107'
    memory "7 GB"

    input:
    tuple val(library), path(metrics)

    output:
    tuple val(library), path("${library}.txt")

    """
    add-metrics.py $metrics > ${library}.txt
    """

}


process plot_qc_metrics {

    publishDir "${params.results}/ataqv/single-nucleus"
    time '10h'
    tag "${library}"
    container 'library://porchard/default/dropkick:20220225'
    memory { 10.GB * task.attempt }
    maxRetries 1
    errorStrategy 'retry'
    cpus 1

    input:
    tuple val(library), path(metrics)

    output:
    tuple val(library), path("*.png")
    path("*.tsv")

    """
    plot-qc-metrics-atac.py --prefix ${library}. $metrics
    """

}


process ATAQV_BULK {
    errorStrategy 'retry'
    maxRetries 1
    memory '8 GB'
    time '8h'
    tag "${sample}"
    container 'docker://porchard/ataqv:1.5.0'
    cpus 2
    publishDir "${params.results}/atac_qc/bulk", mode: 'link'

    input:
    tuple val(sample), val(location), val(cluster_res), val(df_pk)
    
    output:
    tuple val(sample), path("${sample}.bulk.ataqv.json.gz"), emit: json
    tuple val(sample), path("${sample}.bulk.ataqv.txt"), emit: txt

    script:
    """
    ataqv --name ${sample} --ignore-read-groups --metrics-file ${sample}.bulk.ataqv.json.gz --tss-file ${get_tss(params.genome)} ${get_organism(params.genome)} ${location}/atac_possorted_bam.bam > ${sample}.bulk.ataqv.txt
    """
}

process ATAQV_BULK_VIEWER {
    memory '4 GB'
    time '4h'
    tag "${sample}"
    container 'docker://porchard/ataqv:1.5.0'
    cpus 1
    publishDir "${params.results}/atac_qc/bulk_viewer", mode: 'link'

    input:
    tuple val(sample), path(json_metrics)
    
    output:
    tuple val(sample), path("${sample}.html"), emit: viewer

    script:
    """
    # mkarv syntax: mkarv [options] OUTPUT_DIRECTORY INPUT_METRICS.json
    mkarv --force ${sample}.html ${json_metrics}
    """
}

process JOINT_QC {
    tag "${sample_id}"
    memory '32 GB'
    publishDir "${params.results}/joint_qc", mode: 'copy'
    container 'docker://ncarrut/singlecell_qc:second'

    input:
        tuple val(sample_id), val(location), val(cluster_res), val(df_pk), path(atac_metrics), path(intron_counts), path(cellbender_h5), path(rna_metrics), path(knee), path(pass_qc)

    output:
        path "${sample_id}_qcPlot.png", emit: qc_plot
        path "${sample_id}_upsetPlot.png", emit: upset_plot
        path "${sample_id}_metrics.txt", emit: metrics

    script:
    """
    python ${baseDir}/bin/joint_qc.py \\
        --sample ${sample_id} \\
        --ATAC_metrics ${atac_metrics} \\
        --INTRON_COUNTS ${intron_counts} \\
        --CELLBENDER ${cellbender_h5} \\
        --RNA_METRICS ${rna_metrics} \\
        --knee ${knee} \\
        --passQC ${pass_qc} \\
        --RNA_BARCODE_WHITELIST ${params.rna_barcode_whitelist} \\
        --ATAC_BARCODE_WHITELIST ${params.atac_barcode_whitelist} \\
        --qcPlot ${sample_id}_qcPlot.png \\
        --upsetPlot ${sample_id}_upsetPlot.png \\
        --outmetrics ${sample_id}_metrics.txt
    """
}

workflow {
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def sample = row.sample
            def location = row.location
            def cluster_res = row.cluster_res ? row.cluster_res.toFloat() : params.cluster_res
            def df_pk = row.df_pk ? row.df_pk.toFloat() : params.df_pk
            return tuple(sample, location, cluster_res, df_pk)
        }.set { samples_ch }

    samples_ch.view { sample, location, cluster_res, df_pk ->
        "Sample: ${sample} | Location: ${location} | Cluster_res: ${cluster_res} | DoubletFinder_PK: ${df_pk}"
    }
    
    splitter_out = SPLITTER(samples_ch)
    intron_counter_out = INTRONCOUNTER(samples_ch)
    rna_metrics_out = QC(samples_ch.join(splitter_out.GEX))
    qc_out = PLOTQC(rna_metrics_out)
    rankplot_out = INTERACTIVEBARCODERANKPLOT(splitter_out.GEX)
    cellbender_out = CELLBENDER(splitter_out.GEX)
    emptyDrops_out = EMPTYDROPS(splitter_out.GEX.join(cellbender_out.metrics))
    
    // ATAC QC processing with 10X cellranger alignment
    bigwig_out = BIGWIG(samples_ch)
    tss_plot_out = PLOT_SIGNAL_AT_TSS(bigwig_out)
    atac_single_nucleus = ATAQV_SINGLE_NUCLEUS(samples_ch)
    atac_processed = atac_single_nucleus.metrics | add_qc_metrics
    atac_processed | plot_qc_metrics
    atac_bulk = ATAQV_BULK(samples_ch)
    atac_viewer = ATAQV_BULK_VIEWER(atac_bulk.json)

    JOINT_QC(samples_ch.join(atac_processed).join(intron_counter_out.counts).join(cellbender_out.h5_fpr05).join(rna_metrics_out).join(emptyDrops_out.knee_pass))
}
