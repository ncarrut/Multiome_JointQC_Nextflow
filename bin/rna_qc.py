#!/usr/bin/env python3
# coding: utf-8

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import upsetplot

import argparse

from helper_joint_qc import *
from rna_core import compute_rna_metrics

from logging_config import setup_logging
import logging

parser = argparse.ArgumentParser("Plot RNA-only QC metrics per sample (scRNA or snRNA)")
parser.add_argument("--sample", help="Sample ID.", type=str)
parser.add_argument("--RNA_results_dir", help="Path to RNA results directory.", type=str)
parser.add_argument("--assay", help="Assay type.", type=str, choices=["scRNA", "snRNA"])
parser.add_argument("--qcPlot", help="Path to save qcPlot plots.", type=str)
parser.add_argument("--upsetPlot", help="Path to save upset plots.", type=str)
parser.add_argument("--outmetrics", help="Path to save all metrics results.", type=str)
parser.add_argument("--outlogs", help="Path to save log messages.", type=str)

args = parser.parse_args()

# save logs
setup_logging(log_file=args.outlogs, level=logging.DEBUG)
logger = logging.getLogger(__name__)

# ---inputs---
donor = args.sample
RNA_results_dir = str(args.RNA_results_dir)
assay = args.assay
unit_label = "cells" if assay == "scRNA" else "nuclei"

# ---process RNA-side inputs and compute RNA thresholds/filters---
metrics, thresholds, knee_plot_info = compute_rna_metrics(donor, RNA_results_dir, assay=assay)

THRESHOLD_RNA_MIN_UMI = thresholds["rna_min_umi"]
THRESHOLD_FRACTION_CB_REMOVED = thresholds["fraction_cb_removed"]
THRESHOLD_RNA_MAX_MITO = thresholds["rna_max_mito"]

knee = knee_plot_info["knee"]
inflection = knee_plot_info["inflection"]
end_cliff = knee_plot_info["end_cliff"]
plateau = knee_plot_info["plateau"]
n_peaks_knee_plot = knee_plot_info["n_peaks_knee_plot"]

### get cells that passed all thresholds; all RNA-side filters were already assigned by compute_rna_metrics()
metrics['pass_all_filters'] = metrics.filter(like='filter_').all(axis=1)

log_thresholds(thresholds)

##########
# List of pass-QC barcodes
pass_qc_barcodes = list(sorted(metrics[metrics.pass_all_filters].barcode.to_list()))

# Plot QC metrics
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(3*4, 2*4))

ax = axs[0, 0]
barcode_rank_plot(metrics, ax)
ax.axhline(knee, color='red', ls='--', label='knee={:,}'.format(knee))
ax.axhline(inflection, color='green', ls='--', label='inflection={:,}'.format(inflection))
ax.axhline(end_cliff, color='blue', ls='--', label='end_cliff={:,}'.format(end_cliff))
ax.axhline(plateau, color='orange', ls='--', label='plateau={:,}'.format(plateau))
ax.set_title('Inferred n knees = {:,}'.format(n_peaks_knee_plot))
ax.legend()

ax = axs[0, 1]
rna_umis_vs_rna_mito_plot(metrics, ax)
ax.axhline(THRESHOLD_RNA_MAX_MITO/100, color='blue', ls='--', label='THRESHOLD_RNA_MAX_MITO = {:,}'.format(THRESHOLD_RNA_MAX_MITO))
ax.axvline(THRESHOLD_RNA_MIN_UMI, color='red', ls='--')
ax.legend()

ax = axs[0, 2]
cellbender_fraction_removed(metrics, ax)
ax.axhline(THRESHOLD_FRACTION_CB_REMOVED, color='blue', ls='--')

ax = axs[1, 0]
sns.histplot(x='pct_cellbender_removed', data=metrics[(metrics.pct_cellbender_removed > 5) &
                                                      (metrics.pct_cellbender_removed < 50) &
                                                      (np.isnan(metrics.pct_cellbender_removed) == False)], ax=ax)
ax.axvline(THRESHOLD_FRACTION_CB_REMOVED*100, color='blue', ls='--', label='%ambient removed threshold Multi-otsu= {:,}'.format(round(THRESHOLD_FRACTION_CB_REMOVED*100, 2)))
ax.legend()
ax.set_xlabel('5% < % ambient removed < 50%')

ax = axs[1, 1]
cellbender_cell_probabilities(metrics, ax)

ax = axs[1, 2]
if assay == "scRNA":
    THRESHOLD_HELM = thresholds["helm"]
    # same population used by get_helm_threshold() in helper_joint_qc.py
    helm_used_for_threshold = metrics[(metrics.filter_rna_emptyDrops == True) &
                                       (metrics.filter_rna_min_umi == True) &
                                       (metrics.filter_pct_cellbender_removed == True) &
                                       (metrics.filter_rna_max_mito == True) &
                                       (metrics.rna_helm_metric > 0) &
                                       np.isfinite(metrics.rna_log_helm_metric)]
    sns.histplot(x='rna_log_helm_metric', data=helm_used_for_threshold, ax=ax)
    ax.axvline(THRESHOLD_HELM, color='red', ls='--', label='HELM threshold = {:,.2f}'.format(THRESHOLD_HELM))
    ax.set_xlabel('log(HELM) = log(mito frac. x (1 - exon/full-gene ratio))')
    ax.set_title('HELM (intron x mito) metric')
    ax.legend()
else:
    THRESHOLD_EXON_GENE_BODY_RATIO = thresholds["exon_gene_body_ratio"]
    rna_umis_vs_exon_to_full_gene_body_ratio(metrics, ax)
    ax.axhline(THRESHOLD_EXON_GENE_BODY_RATIO, color='red', ls='--', label='exon/full ratio. Multi-otsu = {:,}'.format(round(THRESHOLD_EXON_GENE_BODY_RATIO, 2)))
    ax.legend()
    ax.axvline(THRESHOLD_RNA_MIN_UMI, color='red', ls='--')
    ax.set_xlim(left=0.8*THRESHOLD_RNA_MIN_UMI)

fig.suptitle('{:,} pass QC {} '.format(len(pass_qc_barcodes), unit_label) + donor)
fig.tight_layout()
fig.savefig(args.qcPlot, bbox_inches='tight', dpi=300)

# Plot the number of barcodes passing each filter
fig, ax = plt.subplots(figsize=(7, 6))
ax.remove()

for_upset = metrics.filter(like='filter_').rename(columns=lambda x: 'pass_' + x)
for_upset = for_upset.groupby(for_upset.columns.to_list()).size()
upsetplot.plot(for_upset, fig=fig, sort_by='cardinality', show_counts=True)
fig.savefig(args.upsetPlot, bbox_inches='tight', dpi=300)

metrics.to_csv(args.outmetrics, index=False)
