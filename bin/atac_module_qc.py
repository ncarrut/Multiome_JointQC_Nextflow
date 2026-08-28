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
import tables
import anndata
from typing import Dict, Optional
import numpy as np
import scipy.sparse as sp
from scipy import io
import glob
import os
import upsetplot
from scipy.io import mmread
import csv

import argparse

# import helper functions
try:
    from helper_joint_qc import *
except ModuleNotFoundError as error:
    if error.name == "helper_joint_qc":
        raise SystemExit(
            "Error: Could not find the 'helper_joint_qc' module. "
            "Run this script from the directory containing 'helper_joint_qc.py', "
            "or add that directory's absolute path to Python's module search path "
            "(sys.path or PYTHONPATH)."
        ) from error
    raise

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

parser = argparse.ArgumentParser("Plot QC metrics per sample for snATAC-seq data.")
parser.add_argument("--sample", help="Sample ID. Example: 14806-AH-10-hg38, i.e. `<library name>-<genome>`.", type=str)
parser.add_argument("--ATAC_results_dir", help="Path to ATAC results directory. Path should have subdirectories with upstream analyses such as `ataqv`.", type=str)
parser.add_argument("--lower", help="Read number threshold, a barcode with less than this number is considered empty. Used for min HQAA detection.", type=float, default=100)
parser.add_argument("--filter_MT_ATAC", help="Whether to filter ATAC nuclei based on %chrMT threshold. Default: True.", type=str2bool, default=True)
parser.add_argument("--filter_max_frac_auto_ATAC", help="Whether to filter ATAC nuclei based on maximum fraction of autosomal reads derived from a single autosome threshold. Default: True.", type=str2bool, default=True)
parser.add_argument("--output", help="Path to save QC plots and metrics.", type=str)


args = parser.parse_args()

# ---inputs---
sample = args.sample
logger.info(f"Sample name: {sample}")
ATAC_results_dir = args.ATAC_results_dir
logger.info(f"Input dir for ATAC: {ATAC_results_dir}")
ATAC_METRICS = ATAC_results_dir+'/ataqv/single-nucleus/'+sample+'.txt'

# --- define output files---
output = args.output
qcPlot = output+"/"+sample+".qcPlot.png"
upsetPlot = output+"/"+sample+".upsetPlot.png"
outmetrics = output+"/"+sample+".outmetrics.csv"
outlogs = output+"/"+sample+".log"
if not os.path.isdir(output):
    raise FileNotFoundError(f"Directory does not exist: {output}.")

# save logs
setup_logging(log_file=outlogs, level=logging.DEBUG)
logger = logging.getLogger(__name__)

# ---upfront thresholds--- 
THRESHOLD_ATAC_MIN_TSS_ENRICHMENT = 2

# ---process inputs---
### ATAC side ###
if not os.path.isfile(ATAC_METRICS):
    raise FileNotFoundError(f"File does not exist: {ATAC_METRICS}.")

atac_metrics = pd.read_csv(ATAC_METRICS, sep='\t', index_col=0).rename_axis(index='barcode')
KEEP_ATAC_METRICS = ['median_fragment_length', 'hqaa', 'max_fraction_reads_from_single_autosome', 'percent_mitochondrial', 'tss_enrichment']
atac_metrics = atac_metrics[KEEP_ATAC_METRICS]
atac_metrics.max_fraction_reads_from_single_autosome = atac_metrics.max_fraction_reads_from_single_autosome.fillna(0)
atac_metrics.median_fragment_length = atac_metrics.median_fragment_length.fillna(0)
atac_metrics.percent_mitochondrial = atac_metrics.percent_mitochondrial.fillna(0)
atac_metrics.tss_enrichment = atac_metrics.tss_enrichment.fillna(0)
atac_metrics['fraction_mitochondrial'] = atac_metrics.percent_mitochondrial / 100

metrics = atac_metrics.rename(columns=lambda x: 'atac_' + x)

# get HQAA threshold
values = np.log10(atac_metrics[atac_metrics.hqaa > args.lower].hqaa).values
values = values.reshape((len(values),1))
thresholds = threshold_multiotsu(image=values, classes=3, nbins=256)
# convert back to linear scale
thresholds = [pow(10, i) for i in thresholds]
THRESHOLD_ATAC_MIN_HQAA = round(thresholds[1])

metrics['filter_atac_min_hqaa'] = metrics.atac_hqaa >= THRESHOLD_ATAC_MIN_HQAA

### get THRESHOLD_ATAC_MAX_MITO
if (args.filter_MT_ATAC == True):
    n_peaks, atac_kde_df = guess_n_classes(metrics, "ATAC")
    THRESHOLD_ATAC_MAX_MITO = get_chrMT_threshold_ATAC(metrics, n_peaks = n_peaks)
    logger.info(f"THRESHOLD_ATAC_MAX_MITO = {THRESHOLD_ATAC_MAX_MITO}")

### get THRESHOLD_ATAC_MAX_FRAC_READS_FROM_SINGLE_AUTOSOME
if (args.filter_max_frac_auto_ATAC == True):
    THRESHOLD_ATAC_MAX_FRAC_READS_FROM_SINGLE_AUTOSOME, n_peaks, kde_df = get_atac_max_autosome_threshold(metrics)
    logger.info(f"THRESHOLD_ATAC_MAX_FRAC_READS_FROM_SINGLE_AUTOSOME = {THRESHOLD_ATAC_MAX_FRAC_READS_FROM_SINGLE_AUTOSOME}")


### get cells that passed all thresholds; those that passed post-CB nUMIs have been identified above
metrics['filter_atac_min_hqaa'] = metrics.atac_hqaa >= THRESHOLD_ATAC_MIN_HQAA
metrics['filter_atac_min_tss_enrichment'] = metrics.atac_tss_enrichment >= THRESHOLD_ATAC_MIN_TSS_ENRICHMENT
if (args.filter_max_frac_auto_ATAC == True):
    metrics['filter_max_fraction_reads_from_single_autosome'] = metrics.atac_max_fraction_reads_from_single_autosome <= THRESHOLD_ATAC_MAX_FRAC_READS_FROM_SINGLE_AUTOSOME/100
if (args.filter_MT_ATAC == True):
    metrics['filter_atac_max_mito'] = metrics.atac_percent_mitochondrial <= THRESHOLD_ATAC_MAX_MITO
metrics['pass_all_filters'] = metrics.filter(like='filter_').all(axis=1)

# to collect all Thresholds here
thresholds = {
        "atac_min_hqaa": THRESHOLD_ATAC_MIN_HQAA,
        "atac_min_tss_enrichment": THRESHOLD_ATAC_MIN_TSS_ENRICHMENT
        }

if (args.filter_MT_ATAC == True):
    thresholds["atac_max_mito"] = (THRESHOLD_ATAC_MAX_MITO)

if (args.filter_max_frac_auto_ATAC == True):
    thresholds["atac_max_frac_reads_from_single_autosome"] = (THRESHOLD_ATAC_MAX_FRAC_READS_FROM_SINGLE_AUTOSOME)

log_thresholds(thresholds)


##########
metrics = metrics.reset_index()
# List of pass-QC barcodes
pass_qc_nuclei = list(sorted(metrics[metrics.pass_all_filters].barcode.to_list()))

# Plot QC metrics
fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(2*4, 2*4))

ax=axs[0, 0]
atac_hqaa_vs_atac_tss_enrichment_plot(metrics, ax, alpha=0.02, s=3)
ax.axvline(THRESHOLD_ATAC_MIN_HQAA, color='red', ls='--', label='Min HQAA = {:,}'.format(THRESHOLD_ATAC_MIN_HQAA))
ax.axhline(THRESHOLD_ATAC_MIN_TSS_ENRICHMENT, color='blue', ls='--', label='Min TSS enrichment = {:,}'.format(THRESHOLD_ATAC_MIN_TSS_ENRICHMENT))
ax.legend()

ax=axs[0, 1]
barcode_rank_plot_atac(metrics, ax, alpha=0.02, s=3)
ax.axhline(THRESHOLD_ATAC_MIN_HQAA, color='red', ls='--', label='Min HQAA')
ax.legend()

ax=axs[1, 0]
atac_hqaa_vs_atac_mt_pct_plot(metrics, ax, alpha=0.02, s=3)
ax.axvline(THRESHOLD_ATAC_MIN_HQAA, color='red', ls='--', label='Min HQAA')
if (args.filter_MT_ATAC == True):
    ax.axhline(THRESHOLD_ATAC_MAX_MITO, color='blue', ls='--', label='Max %chrMT = {:,}'.format(THRESHOLD_ATAC_MAX_MITO))
ax.legend()

ax=axs[1, 1]
sns.scatterplot(x='atac_hqaa', y='atac_max_fraction_reads_from_single_autosome', ax=ax, data=metrics, hue='pass_all_filters', palette={True: 'red', False: 'black'}, edgecolor=None, alpha=0.02, s=3)
ax.set_xscale('log')
ax.set_xlabel('HQAA')
ax.set_ylabel('Max fraction reads from single autosome')
ax.axvline(THRESHOLD_ATAC_MIN_HQAA, color='red', ls='--', label='Min HQAA')
if (args.filter_max_frac_auto_ATAC == True):
    ax.axhline(THRESHOLD_ATAC_MAX_FRAC_READS_FROM_SINGLE_AUTOSOME/100, color='blue', ls='--', label='Max frac. from single autosome = {:,}'.format(THRESHOLD_ATAC_MAX_FRAC_READS_FROM_SINGLE_AUTOSOME/100))
ax.legend()

fig.suptitle('{:,} pass-QC nuclei'.format(len(pass_qc_nuclei)) + " " + sample)
fig.tight_layout()
fig.savefig(qcPlot, bbox_inches='tight', dpi=300)

# Plot the number of cells passing each filter
fig, ax = plt.subplots(figsize=(7, 6))
ax.remove()

for_upset = metrics.filter(like='filter_').rename(columns=lambda x: 'pass_' + x)
for_upset = for_upset.groupby(for_upset.columns.to_list()).size()
upsetplot.plot(for_upset, fig=fig, sort_by='cardinality', show_counts=True)
fig.savefig(upsetPlot, bbox_inches='tight', dpi=300)

metrics.to_csv(outmetrics, index=False) 

