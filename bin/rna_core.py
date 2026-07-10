#!/usr/bin/env python
# coding: utf-8

import csv
import os
import logging

import numpy as np
import pandas as pd
from scipy.io import mmread

from helper_joint_qc import *

logger = logging.getLogger(__name__)

THRESHOLD_CELLBENDER_MIN_CELL_PROBABILITY = 0.99

ASSAYS_WITH_EXON_RATIO_FILTER = ("multiome", "snRNA")


def compute_rna_metrics(sample, RNA_results_dir, assay="multiome"):
    """
    Compute RNA-side QC metrics, thresholds, and filters shared by multiome,
    snRNA-alone, and scRNA-alone QC runs.

    For multiome and snRNA, intron retention is gated by the exon/full-gene-body
    ratio filter. For scRNA, that filter is replaced by HELM (mito fraction x
    intron fraction), since whole-cell preps don't carry the same nuclear
    pre-mRNA contamination signature that the exon-ratio filter targets.

    Parameters
    ----------
    sample : str
        Sample/donor ID, used to locate per-sample files under RNA_results_dir.
    RNA_results_dir : str
        Path to the upstream snRNAseq-NextFlow results directory.
    assay : str
        One of "multiome", "snRNA", "scRNA".

    Returns
    -------
    tuple
        (metrics, thresholds, knee_plot_info)
        - metrics : pd.DataFrame, 'barcode' as a column (not the index)
        - thresholds : dict of computed threshold values
        - knee_plot_info : dict with keys knee, inflection, end_cliff, plateau,
          n_peaks_knee_plot
    """
    donor = sample

    logger.info(f"Sample name: {donor}")
    logger.info(f"Input dir for RNA: {RNA_results_dir}")
    logger.info(f"Assay: {assay}")

    CELLBENDER = RNA_results_dir + 'cellbender/' + donor + '-hg38.cellbender_FPR_0.05.h5'
    RNA_METRICS = RNA_results_dir + 'qc/' + donor + '-hg38.qc.txt'
    GENE_FULL_EXON_OVER_INTRON_COUNTS = RNA_results_dir + 'starsolo/' + donor + '-hg38/' + donor + '-hg38.Solo.out/GeneFull_ExonOverIntron/raw'
    GENE_COUNTS = RNA_results_dir + 'starsolo/' + donor + '-hg38/' + donor + '-hg38.Solo.out/Gene/raw'
    knee_file = RNA_results_dir + 'emptyDrops/' + donor + '-hg38.knee.txt'
    passQC = RNA_results_dir + 'emptyDrops/' + donor + '-hg38.pass.txt'

    ## load metrics df
    adata = anndata_from_h5(CELLBENDER, analyzed_barcodes_only=True)
    rna_metrics = pd.read_csv(RNA_METRICS, sep='\t')
    rna_metrics = rna_metrics[rna_metrics.barcode != '-']

    ## Calculate ratio of exonic vs full gene body reads
    # exons only
    gene_mat = mmread(os.path.join(GENE_COUNTS, 'matrix.mtx'))
    gene_umis_per_barcode = gene_mat.sum(axis=0).tolist()[0]

    # includes introns
    gene_full_mat = mmread(os.path.join(GENE_FULL_EXON_OVER_INTRON_COUNTS, 'matrix.mtx'))
    gene_full_umis_per_barcode = gene_full_mat.sum(axis=0).tolist()[0]

    barcodes = pd.read_csv(os.path.join(GENE_COUNTS, 'barcodes.tsv'), header=None)[0]
    assert(all(barcodes == pd.read_csv(os.path.join(GENE_FULL_EXON_OVER_INTRON_COUNTS, 'barcodes.tsv'), header=None)[0]))

    exon_to_full_gene_body_ratio = pd.DataFrame({'barcode': barcodes, 'gene': gene_umis_per_barcode, 'gene_full': gene_full_umis_per_barcode})
    exon_to_full_gene_body_ratio['exon_to_full_gene_body_ratio'] = exon_to_full_gene_body_ratio.gene / exon_to_full_gene_body_ratio.gene_full
    rna_metrics = rna_metrics.merge(exon_to_full_gene_body_ratio)
    metrics = rna_metrics.set_index('barcode').rename(columns=lambda x: 'rna_' + x)

    ## cellbender-related statistics
    metrics = metrics.reset_index()
    cell_probability = cellbender_anndata_to_cell_probability(adata)
    post_cellbender_umis = umi_count_after_decontamination(adata)

    metrics['cell_probability'] = metrics.barcode.map(lambda x: cell_probability[x] if x in cell_probability else np.nan)
    metrics['post_cellbender_umis'] = metrics.barcode.map(lambda x: post_cellbender_umis[x] if x in post_cellbender_umis else np.nan)
    metrics['fraction_cellbender_removed'] = (metrics.rna_umis - metrics.post_cellbender_umis) / metrics.rna_umis
    metrics['rna_percent_mitochondrial'] = metrics.rna_fraction_mitochondrial * 100
    metrics['pct_cellbender_removed'] = metrics.fraction_cellbender_removed * 100
    metrics['filter_cellbender_cell_probability'] = metrics.cell_probability >= THRESHOLD_CELLBENDER_MIN_CELL_PROBABILITY

    ### get bc that passed emptydrops analysis
    bc = pd.read_csv(passQC, header=0, delim_whitespace="\t")
    metrics['filter_rna_emptyDrops'] = metrics['barcode'].isin(bc.barcode)

    ### load metrics on knee plot
    with open(knee_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        next(reader, None)
        for row in reader:
            knee = round(float(row[0]))
            inflection = round(float(row[1]))
            inflection_rank = round(float(row[2]))
            knee_rank = round(float(row[3]))
            end_cliff = round(float(row[4]))
            end_cliff_rank = round(float(row[5]))
            plateau = round(float(row[6]))

    ### get bc that passed threshold UMIs obtained using Multi-Otsu
    # try to infer UMI threshold
    MAX_EXPECTED_NUMBER_NUCLEI = round_up(len(metrics[metrics.rna_umis >= inflection]), 3)
    LOWERBOUNDS = np.concatenate(([1, 5], np.arange(10, 251, 10), [300, 350, 400, 450, 500]))

    for i in LOWERBOUNDS:
        UMI_THRESHOLD = estimate_threshold(metrics[(metrics.barcode != '-') & (metrics.rna_umis >= i)].rna_umis.astype(int))
        NUMBER_MEETING_UMI_THRESHOLD = (metrics.rna_umis >= UMI_THRESHOLD).sum()
        #allow 1% wiggle room
        if (NUMBER_MEETING_UMI_THRESHOLD*101/100 <= MAX_EXPECTED_NUMBER_NUCLEI) or (NUMBER_MEETING_UMI_THRESHOLD*99/100 <= MAX_EXPECTED_NUMBER_NUCLEI):
            break

    if (NUMBER_MEETING_UMI_THRESHOLD*101/100 > MAX_EXPECTED_NUMBER_NUCLEI) and (NUMBER_MEETING_UMI_THRESHOLD*99/100 > MAX_EXPECTED_NUMBER_NUCLEI):
        # just fall back to 500
        UMI_THRESHOLD = 500
        NUMBER_MEETING_UMI_THRESHOLD = (metrics.rna_umis>=UMI_THRESHOLD).sum()

    THRESHOLD_RNA_MIN_UMI = UMI_THRESHOLD
    metrics['filter_rna_min_umi'] = metrics.rna_umis >= THRESHOLD_RNA_MIN_UMI

    # get %ambient vs post CB UMI thresholding
    peaks_cb, n_peaks_cb, cb_kde_df = guess_n_classes_cellbender(metrics)
    THRESHOLD_FRACTION_CB_REMOVED, THRESHOLD_POST_CB_UMIS = get_cellbender_thresholds(metrics, peaks_cb, n_peaks_cb, cb_kde_df)

    metrics['filter_pct_cellbender_removed'] = metrics.pct_cellbender_removed <= THRESHOLD_FRACTION_CB_REMOVED*100

    ### get THRESHOLD_EXON_GENE_BODY_RATIO (multiome / snRNA only; scRNA uses HELM instead, computed below)
    THRESHOLD_EXON_GENE_BODY_RATIO = None
    if assay in ASSAYS_WITH_EXON_RATIO_FILTER:
        x = np.log10(metrics[(metrics.rna_exon_to_full_gene_body_ratio>0)&
                             (metrics.filter_rna_min_umi ==True)&
                             (metrics.filter_pct_cellbender_removed ==True)].rna_umis)
        y = metrics[(metrics.rna_exon_to_full_gene_body_ratio>0)&
                    (metrics.filter_rna_min_umi ==True)&
                    (metrics.filter_pct_cellbender_removed ==True)].rna_exon_to_full_gene_body_ratio

        THRESHOLD_EXON_GENE_BODY_RATIO = get_exon_fullgene_ratio(x, y)

        if THRESHOLD_EXON_GENE_BODY_RATIO >= 0.95:
            data = metrics[(metrics.rna_exon_to_full_gene_body_ratio>0)&
                      (metrics.rna_exon_to_full_gene_body_ratio<1.0)&
                      (metrics.filter_rna_min_umi ==True)&
                      (metrics.filter_pct_cellbender_removed ==True)].rna_exon_to_full_gene_body_ratio.astype(float).values
            THRESHOLD_EXON_GENE_BODY_RATIO = threshold_multiotsu(data, classes=3)[1]

    ### get THRESHOLD_RNA_MAX_MITO
    n_peaks, rna_kde_df = guess_n_classes(metrics, "RNA")
    THRESHOLD_RNA_MAX_MITO = get_chrMT_threshold_RNA(metrics, n_peaks = n_peaks)

    ### assign filters that depend on the thresholds computed above
    metrics['filter_rna_max_mito'] = metrics.rna_percent_mitochondrial <= THRESHOLD_RNA_MAX_MITO

    THRESHOLD_HELM = None
    if assay in ASSAYS_WITH_EXON_RATIO_FILTER:
        metrics['filter_rna_exon_to_full_gene_body_ratio'] = metrics.rna_exon_to_full_gene_body_ratio <= THRESHOLD_EXON_GENE_BODY_RATIO
    else:
        # scRNA: HELM = mito fraction x intron fraction, replaces the exon-ratio filter
        metrics['rna_helm_metric'] = metrics.rna_fraction_mitochondrial * (1 - metrics.rna_exon_to_full_gene_body_ratio)
        with np.errstate(divide='ignore', invalid='ignore'):
            metrics['rna_log_helm_metric'] = np.log(metrics['rna_helm_metric'])
        THRESHOLD_HELM, n_peaks_helm = get_helm_threshold(metrics)
        metrics['filter_helm'] = metrics.rna_log_helm_metric >= THRESHOLD_HELM

    ##############################

    ####### knee plot analysis
    df_ranked, df_interpolated, n_peaks_knee_plot, final_peak_indices = analyze_knee_plot(metrics, knee, knee_rank, end_cliff, end_cliff_rank, inflection_rank)
    ##############################

    thresholds = {
        "rna_min_umi": THRESHOLD_RNA_MIN_UMI,
        "fraction_cb_removed": THRESHOLD_FRACTION_CB_REMOVED,
        "rna_max_mito": THRESHOLD_RNA_MAX_MITO,
    }
    if assay in ASSAYS_WITH_EXON_RATIO_FILTER:
        thresholds["exon_gene_body_ratio"] = THRESHOLD_EXON_GENE_BODY_RATIO
    else:
        thresholds["helm"] = THRESHOLD_HELM

    knee_plot_info = {
        "knee": knee,
        "inflection": inflection,
        "end_cliff": end_cliff,
        "plateau": plateau,
        "n_peaks_knee_plot": n_peaks_knee_plot,
    }

    return metrics, thresholds, knee_plot_info
