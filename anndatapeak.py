import pdb,sys,os
import anndata
import scanpy as sc
#from File import *
import pandas as pd
import argparse
import matplotlib
matplotlib.use('Agg')

fn = 'C:/Users/Andre/Downloads/PPUCMUCPT_out10224/PPUCMUCPTfiltered_feature_bc_matrix14.E.h5ad'
#fn = 'C:/Users/Andre/Downloads/example_out1229/example.E.h5ad'

adata_orig = sc.read_h5ad(fn)
adata_orig.obs_names_make_unique()
adata_orig.var_names_make_unique()
#print(adata_orig)
#var_names_make_unique()
adata_orig.obs['leiden'] = adata_orig.obs['leiden'].astype(int)
print(adata_orig.obs.head())
print(adata_orig.obs.dtypes)
adata_orig.write_h5ad('C:/Users/Andre/Downloads/PPUCMUCPT_out10224/PPUCMUCPTfiltered_feature_bc_matrix14_2.E.h5ad',compression=9)