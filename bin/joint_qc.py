#!/usr/bin/env python3
# coding: utf-8

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

import os

current_dir = os.getcwd()
os.environ['MPLCONFIGDIR'] = current_dir

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

parser = argparse.ArgumentParser("Plot QC metrics per sample")
parser.add_argument("--sample", help="Donor ID.", type=str)
parser.add_argument("--ATAC_metrics", help="Path to ATAC ataqv metrics file.", type=str)
parser.add_argument("--INTRON_COUNTS", help="Path to intron counter output file.", type=str)
parser.add_argument("--CELLBENDER", help="Path to CellBender FPR 0.05 h5 file.", type=str)
parser.add_argument("--RNA_METRICS", help="Path to RNA QC metrics file.", type=str)
parser.add_argument("--knee", help="Path to emptyDrops knee file.", type=str)
parser.add_argument("--passQC", help="Path to emptyDrops pass file.", type=str)
parser.add_argument("--filter_MT_ATAC", help="Whether to filter ATAC nuclei based on %chrMT threshold. Default: True.", type=str2bool, default=True)
parser.add_argument("--qcPlot", help="Path to save qcPlot plots.", type=str)
parser.add_argument("--upsetPlot", help="Path to save upset plots.", type=str)
parser.add_argument("--outmetrics", help="Path to save all metrics results.", type=str)
parser.add_argument("--manual", action="store_true", help="Manual mode: read all args from .command.sh in the current working directory.")


args = parser.parse_args()

#args.manual = True
# ---inputs---
if args.manual:
    import re
    with open('.command.sh') as _f:
        _cmd = _f.read()
    def _get_arg(name):
        m = re.search(rf'--{name}\s+(\S+)', _cmd)
        return m.group(1) if m else None
    sample = _get_arg('sample')
    ATAC_METRICS = _get_arg('ATAC_metrics')
    INTRON_COUNTS = _get_arg('INTRON_COUNTS')
    CELLBENDER = _get_arg('CELLBENDER')
    RNA_METRICS = _get_arg('RNA_METRICS')
    knee = _get_arg('knee')
    passQC = _get_arg('passQC')
    args.qcPlot = _get_arg('qcPlot')
    args.upsetPlot = _get_arg('upsetPlot')
    args.outmetrics = _get_arg('outmetrics')
else:
    sample = args.sample
    ATAC_METRICS = args.ATAC_metrics
    INTRON_COUNTS = args.INTRON_COUNTS
    CELLBENDER = args.CELLBENDER
    RNA_METRICS = args.RNA_METRICS
    knee = args.knee
    passQC = args.passQC

# --- define output files---
qcPlot = args.qcPlot
upsetPlot = args.upsetPlot
outmetrics = args.outmetrics
outlogs = os.path.join(os.path.dirname(outmetrics) or '.', sample + '.log')

# save logs
setup_logging(log_file=outlogs, level=logging.DEBUG)
logger = logging.getLogger(__name__)

logger.info(f"Sample name: {sample}")

# sanity check input files
for _label, _path in [("CELLBENDER", CELLBENDER), ("RNA_METRICS", RNA_METRICS),
                       ("ATAC_METRICS", ATAC_METRICS), ("INTRON_COUNTS", INTRON_COUNTS),
                       ("knee", knee), ("passQC", passQC)]:
    if not os.path.isfile(_path):
        raise FileNotFoundError(f"File does not exist: {_path}.")

# ---upfront thresholds---
THRESHOLD_CELLBENDER_MIN_CELL_PROBABILITY = 0.99
THRESHOLD_ATAC_MIN_TSS_ENRICHMENT = 2

# ---process inputs---
## Note: unlike the custom snRNAseq-NextFlow/snATACseq-NextFlow pipelines, 10x
## cellranger-arc already emits RNA and ATAC data under a shared barcode
## namespace, so no ATAC->RNA barcode whitelist translation is needed here.

## load metrics df
adata = anndata_from_h5(CELLBENDER, analyzed_barcodes_only=True)
rna_metrics = pd.read_csv(RNA_METRICS, sep='\t')
rna_metrics = rna_metrics[rna_metrics.barcode!='-']

## Calculate ratio of exonic vs full gene body reads from intron counter output
intron_counts = pd.read_csv(INTRON_COUNTS, sep=r'\s+', header=None, names=['count', 'barcode', 'type'],
                             engine='python', skipinitialspace=True)
intron_counts = intron_counts[intron_counts.type.isin(['E', 'I'])]
intron_pivot = intron_counts.pivot_table(index='barcode', columns='type', values='count',
                                         aggfunc='sum', fill_value=0).reset_index()
intron_pivot = intron_pivot.set_index('barcode')
intron_pivot.columns.name = None

exon_ratio_dict = (intron_pivot['E'] / (intron_pivot['E'] + intron_pivot['I'])).to_dict()
rna_metrics['exon_to_full_gene_body_ratio'] = rna_metrics.barcode.map(exon_ratio_dict)
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
KNEE_FILE = knee
with open(KNEE_FILE, 'r') as file:
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
    UMI_THRESHOLD = estimate_threshold(metrics[(metrics.barcode!='-') & (metrics.rna_umis>=i)].rna_umis.astype(int))
    NUMBER_MEETING_UMI_THRESHOLD = (metrics.rna_umis>=UMI_THRESHOLD).sum()
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


### get THRESHOLD_EXON_GENE_BODY_RATIO
x = np.log10(metrics[(metrics.rna_exon_to_full_gene_body_ratio>0)&
                     (metrics.filter_rna_min_umi ==True)].rna_umis)
y = metrics[(metrics.rna_exon_to_full_gene_body_ratio>0)&
            (metrics.filter_rna_min_umi ==True)].rna_exon_to_full_gene_body_ratio

THRESHOLD_EXON_GENE_BODY_RATIO = get_exon_fullgene_ratio(x, y)

if THRESHOLD_EXON_GENE_BODY_RATIO >= 0.95:
    data = metrics[(metrics.rna_exon_to_full_gene_body_ratio>0)&
              (metrics.rna_exon_to_full_gene_body_ratio<1.0)].rna_exon_to_full_gene_body_ratio.astype(float).values
    THRESHOLD_EXON_GENE_BODY_RATIO = threshold_multiotsu(data, classes=3)[1]

### get THRESHOLD_RNA_MAX_MITO
n_peaks, rna_kde_df = guess_n_classes(metrics, "RNA")
THRESHOLD_RNA_MAX_MITO = get_chrMT_threshold_RNA(metrics, n_peaks = n_peaks)
##############################

####### knee plot analysis
from scipy.interpolate import interp1d
from scipy.signal import find_peaks, savgol_filter

df_ranked, df_interpolated, n_peaks_knee_plot, final_peak_indices = analyze_knee_plot(metrics, knee, knee_rank, end_cliff, end_cliff_rank, inflection_rank)
##############################

### ATAC side ###
atac_metrics = pd.read_csv(ATAC_METRICS, sep='\t', index_col=0).rename_axis(index='barcode')
KEEP_ATAC_METRICS = ['median_fragment_length', 'hqaa', 'max_fraction_reads_from_single_autosome', 'percent_mitochondrial', 'tss_enrichment']
atac_metrics = atac_metrics[KEEP_ATAC_METRICS]
atac_metrics.max_fraction_reads_from_single_autosome = atac_metrics.max_fraction_reads_from_single_autosome.fillna(0)
atac_metrics.median_fragment_length = atac_metrics.median_fragment_length.fillna(0)
atac_metrics.percent_mitochondrial = atac_metrics.percent_mitochondrial.fillna(0)
atac_metrics.tss_enrichment = atac_metrics.tss_enrichment.fillna(0)
atac_metrics['fraction_mitochondrial'] = atac_metrics.percent_mitochondrial / 100

# cellranger-arc's RNA and ATAC barcodes already share one namespace, so no
# atac_to_rna barcode translation is applied here (unlike main's joint_qc.py,
# which is fed by separate custom RNA/ATAC pipelines with distinct whitelists).
metrics = metrics.set_index('barcode').rename(columns=lambda x: '' + x).join(atac_metrics.rename(columns=lambda x: 'atac_' + x))

# get HQAA threshold
values = np.log10(atac_metrics[(atac_metrics.tss_enrichment > 2)].hqaa).values
values = values.reshape((len(values),1))
thresholds = threshold_multiotsu(image=values, classes=2, nbins=256)
# convert back to linear scale
thresholds = [pow(10, i) for i in thresholds]
lower_thres = round(thresholds[0])
lower_thres = max(lower_thres, 100)
values = np.log10(atac_metrics[(atac_metrics.hqaa > lower_thres)].hqaa).values
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



### get cells that passed all thresholds; those that passed post-CB nUMIs have been identified above
metrics['filter_cellbender_cell_probability'] = metrics.cell_probability >= THRESHOLD_CELLBENDER_MIN_CELL_PROBABILITY
metrics['filter_rna_max_mito'] = metrics.rna_percent_mitochondrial <= THRESHOLD_RNA_MAX_MITO
metrics['filter_rna_exon_to_full_gene_body_ratio'] = metrics.rna_exon_to_full_gene_body_ratio <= THRESHOLD_EXON_GENE_BODY_RATIO
metrics['filter_atac_min_hqaa'] = metrics.atac_hqaa >= THRESHOLD_ATAC_MIN_HQAA
metrics['filter_atac_min_tss_enrichment'] = metrics.atac_tss_enrichment >= THRESHOLD_ATAC_MIN_TSS_ENRICHMENT
if (args.filter_MT_ATAC == True):
    metrics['filter_atac_max_mito'] = metrics.atac_percent_mitochondrial <= THRESHOLD_ATAC_MAX_MITO
metrics['pass_all_filters'] = metrics.filter(like='filter_').all(axis=1)

# to collect all Thresholds here
if (args.filter_MT_ATAC == True):
    thresholds = {
        "rna_min_umi": THRESHOLD_RNA_MIN_UMI,
        "fraction_cb_removed": THRESHOLD_FRACTION_CB_REMOVED,
        "rna_max_mito": THRESHOLD_RNA_MAX_MITO,
        "exon_gene_body_ratio": THRESHOLD_EXON_GENE_BODY_RATIO,
        "atac_min_hqaa": THRESHOLD_ATAC_MIN_HQAA,
        "atac_min_tss_enrichment": THRESHOLD_ATAC_MIN_TSS_ENRICHMENT,
        "atac_max_mito": THRESHOLD_ATAC_MAX_MITO,
        }
else:
    thresholds = {
        "rna_min_umi": THRESHOLD_RNA_MIN_UMI,
        "fraction_cb_removed": THRESHOLD_FRACTION_CB_REMOVED,
        "rna_max_mito": THRESHOLD_RNA_MAX_MITO,
        "exon_gene_body_ratio": THRESHOLD_EXON_GENE_BODY_RATIO,
        "atac_min_hqaa": THRESHOLD_ATAC_MIN_HQAA,
        "atac_min_tss_enrichment": THRESHOLD_ATAC_MIN_TSS_ENRICHMENT
        }

log_thresholds(thresholds)

knee_plot_info = {
    "knee": knee,
    "inflection": inflection,
    "end_cliff": end_cliff,
    "plateau": plateau,
    "n_peaks_knee_plot": n_peaks_knee_plot,
}
log_thresholds(knee_plot_info)


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
ax.axhline(THRESHOLD_RNA_MAX_MITO/100, color='blue', ls='--', label='max RNA %chrMT = {:,}'.format(THRESHOLD_RNA_MAX_MITO))
ax.axvline(THRESHOLD_RNA_MIN_UMI, color='red', ls='--')
ax.legend(loc="upper right")

ax = axs[0, 2]
cellbender_fraction_removed(metrics, ax)
ax.axhline(THRESHOLD_FRACTION_CB_REMOVED, color='blue', ls='--', label='%ambient removed threshold Multi-otsu= {:,}'.format(THRESHOLD_FRACTION_CB_REMOVED))

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
ax.axvline(THRESHOLD_ATAC_MIN_HQAA, color='red', ls='--', label='min ATAC HQAA = {:,}'.format(THRESHOLD_ATAC_MIN_HQAA))
ax.axhline(THRESHOLD_ATAC_MIN_TSS_ENRICHMENT, color='red', ls='--')
ax.legend()

ax = axs[2, 2]
atac_hqaa_vs_atac_mt_pct_plot(metrics, ax, alpha=0.02)
ax.axvline(THRESHOLD_ATAC_MIN_HQAA, color='red', ls='--')
if (args.filter_MT_ATAC == True):
    ax.axhline(THRESHOLD_ATAC_MAX_MITO, color='green', ls='--', label='max ATAC %chrMT = {:,}'.format(THRESHOLD_ATAC_MAX_MITO))
ax.legend()


fig.suptitle('{:,} pass QC nuclei'.format(len(pass_qc_nuclei)) + " " + sample)
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
