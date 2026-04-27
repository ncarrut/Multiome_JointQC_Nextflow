#!/usr/bin/env python3

import scanpy as sc
import argparse
import scipy.io as sio
import pandas as pd
import os
import gzip
import shutil


def write_10x_mtx(adata, out_dir, layer=None):
    os.makedirs(out_dir, exist_ok=True)
    # Get matrix (transpose to genes x cells)
    mat = adata.layers[layer].T if layer else adata.X.T
    sio.mmwrite(os.path.join(out_dir, "matrix.mtx"), mat)
    pd.DataFrame(adata.obs_names).to_csv(os.path.join(out_dir, "barcodes.tsv"),sep="\t", header=False, index=False)
    adata.var[["gene_ids", "feature_types"]].reset_index().to_csv(os.path.join(out_dir, "features.tsv"),sep="\t", header=False, index=False)
    #cellbender wants features to be called genes
    adata.var[["gene_ids", "feature_types"]].reset_index().to_csv(os.path.join(out_dir, "genes.tsv"),sep="\t", header=False, index=False)



parser = argparse.ArgumentParser("Split 10 raw h5 file into GEX and ATAC components")
parser.add_argument("--h5", help="h5 location", type=str)
parser.add_argument("--sample", help="Donor ID.", type=str)
args = parser.parse_args()

h5 = args.h5 
sample = args.sample

adata = sc.read_10x_h5(h5, gex_only=False)

adata.var_names_make_unique()
rna = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()
atac = adata[:, adata.var["feature_types"] == "Peaks"].copy()

write_10x_mtx(rna, out_dir=sample + '_GEX')
write_10x_mtx(atac, out_dir=sample + '_ATAC')

