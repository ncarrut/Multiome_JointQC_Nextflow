#!/usr/bin/env python3
# coding: utf-8

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import numpy as np
import os
import upsetplot

import argparse

from helper_joint_qc import *
from rna_core import compute_rna_metrics

from logging_config import setup_logging
import logging

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

parser = argparse.ArgumentParser("Plot QC metrics per sample")
parser.add_argument("--sample", help="Sample ID.", type=str)
parser.add_argument("--RNA_results_dir", help="Path to RNA results directory.", type=str)
parser.add_argument("--ATAC_results_dir", help="Path to ATAC results directory.", type=str)
parser.add_argument("--RNA_BARCODE_WHITELIST", help="Path to RNA barcode whitelist.", type=str)
parser.add_argument("--ATAC_BARCODE_WHITELIST", help="Path to ATAC barcode whitelist.", type=str)
parser.add_argument("--filter_MT_ATAC", help="Whether to filter ATAC nuclei based on %chrMT threshold. Default: False.", type=str2bool, default=False)
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
ATAC_results_dir = args.ATAC_results_dir
logger.info(f"Input dir for ATAC: {ATAC_results_dir}")
RNA_BARCODE_WHITELIST = args.RNA_BARCODE_WHITELIST
ATAC_BARCODE_WHITELIST = args.ATAC_BARCODE_WHITELIST

# ---upfront thresholds---
THRESHOLD_ATAC_MIN_TSS_ENRICHMENT = 2

# ---process RNA-side inputs and compute RNA thresholds/filters---
metrics, thresholds, knee_plot_info = compute_rna_metrics(donor, RNA_results_dir, assay="multiome")

THRESHOLD_RNA_MIN_UMI = thresholds["rna_min_umi"]
THRESHOLD_FRACTION_CB_REMOVED = thresholds["fraction_cb_removed"]
THRESHOLD_RNA_MAX_MITO = thresholds["rna_max_mito"]
THRESHOLD_EXON_GENE_BODY_RATIO = thresholds["exon_gene_body_ratio"]

knee = knee_plot_info["knee"]
inflection = knee_plot_info["inflection"]
end_cliff = knee_plot_info["end_cliff"]
plateau = knee_plot_info["plateau"]
n_peaks_knee_plot = knee_plot_info["n_peaks_knee_plot"]

## ATAC --> RNA barcode mappings
rna_barcodes = pd.read_csv(RNA_BARCODE_WHITELIST, header=None)[0].to_list()
atac_barcodes = pd.read_csv(ATAC_BARCODE_WHITELIST, header=None)[0].to_list()
atac_to_rna = dict(zip(atac_barcodes, rna_barcodes))

### ATAC side ###
atac_metrics = pd.read_csv(ATAC_results_dir+'ataqv/single-nucleus/'+donor+'-hg38.txt', sep='\t', index_col=0).rename_axis(index='barcode')
KEEP_ATAC_METRICS = ['median_fragment_length', 'hqaa', 'max_fraction_reads_from_single_autosome', 'percent_mitochondrial', 'tss_enrichment']
atac_metrics = atac_metrics[KEEP_ATAC_METRICS]
atac_metrics.max_fraction_reads_from_single_autosome = atac_metrics.max_fraction_reads_from_single_autosome.fillna(0)
atac_metrics.median_fragment_length = atac_metrics.median_fragment_length.fillna(0)
atac_metrics.percent_mitochondrial = atac_metrics.percent_mitochondrial.fillna(0)
atac_metrics.tss_enrichment = atac_metrics.tss_enrichment.fillna(0)
atac_metrics['fraction_mitochondrial'] = atac_metrics.percent_mitochondrial / 100

atac_metrics.index = atac_metrics.index.map(atac_to_rna)

metrics = metrics.set_index('barcode').rename(columns=lambda x: '' + x).join(atac_metrics.rename(columns=lambda x: 'atac_' + x))

# get HQAA threshold
values = np.log10(atac_metrics[(atac_metrics.tss_enrichment > 2)].hqaa).values
values = values.reshape((len(values),1))
thresholds_multiotsu = threshold_multiotsu(image=values, classes=2, nbins=256)
# convert back to linear scale
thresholds_multiotsu = [pow(10, i) for i in thresholds_multiotsu]
lower_thres = round(thresholds_multiotsu[0])
lower_thres = max(lower_thres, 100)
values = np.log10(atac_metrics[(atac_metrics.hqaa > lower_thres)].hqaa).values
values = values.reshape((len(values),1))
thresholds_multiotsu = threshold_multiotsu(image=values, classes=3, nbins=256)
# convert back to linear scale
thresholds_multiotsu = [pow(10, i) for i in thresholds_multiotsu]
THRESHOLD_ATAC_MIN_HQAA = round(thresholds_multiotsu[1])

metrics['filter_atac_min_hqaa'] = metrics.atac_hqaa >= THRESHOLD_ATAC_MIN_HQAA

### get THRESHOLD_ATAC_MAX_MITO
if (args.filter_MT_ATAC == True):
    n_peaks, atac_kde_df = guess_n_classes(metrics, "ATAC")
    THRESHOLD_ATAC_MAX_MITO = get_chrMT_threshold_ATAC(metrics, n_peaks = n_peaks)



### get cells that passed all thresholds; RNA-side filters were already assigned by compute_rna_metrics()
metrics['filter_atac_min_hqaa'] = metrics.atac_hqaa >= THRESHOLD_ATAC_MIN_HQAA
metrics['filter_atac_min_tss_enrichment'] = metrics.atac_tss_enrichment >= THRESHOLD_ATAC_MIN_TSS_ENRICHMENT
if (args.filter_MT_ATAC == True):
    metrics['filter_atac_max_mito'] = metrics.atac_percent_mitochondrial <= THRESHOLD_ATAC_MAX_MITO
metrics['pass_all_filters'] = metrics.filter(like='filter_').all(axis=1)

# to collect all Thresholds here
thresholds["atac_min_hqaa"] = THRESHOLD_ATAC_MIN_HQAA
thresholds["atac_min_tss_enrichment"] = THRESHOLD_ATAC_MIN_TSS_ENRICHMENT
if (args.filter_MT_ATAC == True):
    thresholds["atac_max_mito"] = THRESHOLD_ATAC_MAX_MITO

log_thresholds(thresholds)


##########
metrics = metrics.reset_index()
# List of pass-QC barcodes
pass_qc_nuclei = list(sorted(metrics[metrics.pass_all_filters].barcode.to_list()))


# Plot QC metrics #to work on plotting
fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(3*4, 3*4))

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
rna_umis_vs_exon_to_full_gene_body_ratio(metrics, ax)
ax.axhline(THRESHOLD_EXON_GENE_BODY_RATIO, color='red', ls='--', label='exon/full ratio. Multi-otsu = {:,}'.format(round(THRESHOLD_EXON_GENE_BODY_RATIO, 2)))
ax.legend()
ax.axvline(THRESHOLD_RNA_MIN_UMI, color='red', ls='--')
ax.set_xlim(left=0.8*THRESHOLD_RNA_MIN_UMI)

ax = axs[2, 0]
rna_umis_vs_atac_hqaa_plot(metrics, ax)
ax.axhline(THRESHOLD_ATAC_MIN_HQAA, color='red', ls='--')
ax.axvline(THRESHOLD_RNA_MIN_UMI, color='red', ls='--')

ax = axs[2, 1]
atac_hqaa_vs_atac_tss_enrichment_plot(metrics, ax, alpha=0.02)
ax.axvline(THRESHOLD_ATAC_MIN_HQAA, color='red', ls='--', label='THRESHOLD_ATAC_MIN_HQAA = {:,}'.format(THRESHOLD_ATAC_MIN_HQAA))
ax.axhline(THRESHOLD_ATAC_MIN_TSS_ENRICHMENT, color='red', ls='--')
ax.legend()

ax = axs[2, 2]
#barcode_rank_plot_atac(metrics, ax, alpha=0.02)
#ax.axhline(THRESHOLD_ATAC_MIN_HQAA, color='red', ls='--')

#ax = axs[2, 3]
atac_hqaa_vs_atac_mt_pct_plot(metrics, ax, alpha=0.02)
ax.axvline(THRESHOLD_ATAC_MIN_HQAA, color='red', ls='--')
if (args.filter_MT_ATAC == True):
    ax.axhline(THRESHOLD_ATAC_MAX_MITO, color='green', ls='--', label='THRESHOLD_ATAC_MAX_MITO = {:,}'.format(THRESHOLD_ATAC_MAX_MITO))
ax.legend()


fig.suptitle('{:,} pass QC nuclei'.format(len(pass_qc_nuclei)) + " " + donor)
fig.tight_layout()
fig.savefig(args.qcPlot, bbox_inches='tight', dpi=300)

# Plot the number of cells passing each filter
fig, ax = plt.subplots(figsize=(7, 6))
ax.remove()

for_upset = metrics.filter(like='filter_').rename(columns=lambda x: 'pass_' + x)
for_upset = for_upset.groupby(for_upset.columns.to_list()).size()
upsetplot.plot(for_upset, fig=fig, sort_by='cardinality', show_counts=True)
fig.savefig(args.upsetPlot, bbox_inches='tight', dpi=300)


metrics.to_csv(args.outmetrics, index=False)
