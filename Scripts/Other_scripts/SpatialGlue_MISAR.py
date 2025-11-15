
#SpatialGlue for RNA

#data integration for Spatial-ATAC-RNA-seq data
#if using jupyter notebook, select python39 kernel

import os
print(os.getcwd())

import sys
print(sys.version)

import os
import torch
import pandas as pd
import seaborn as sns
import scanpy as sc
import squidpy as sq
from SpatialGlue import SpatialGlue
import matplotlib.pyplot as plt
import argparse
import re

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--sampleId", required=True, help="input_image")

args = vars(ap.parse_args())

sampleId = args["sampleId"]


# Environment configuration. SpatialGlue pacakge can be implemented with either CPU or GPU. GPU acceleration is highly recommend for imporoved efficiency.
device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')

# the location of R, which is necessary for mclust algorithm. Please replace the path below with local R installation path
os.environ['R_HOME'] = '/public/home/fqjiang_gdl/anaconda3/envs/python37/lib/R'

# read data
adata_RNA = sc.read_h5ad(sampleId + '-RNA.h5ad')
adata_ATAC = sc.read_h5ad( sampleId + '-ATAC.h5ad')

#if the data haven't adata.obsm['spatial'] spatial location info.
#create adata.obsm['spatial']
#cord_RNA=pd.DataFrame(data=adata_RNA.obs, index=adata_RNA.obs_names, columns=['array_col','array_row'])
#adata_RNA.obsm['spatial'] = cord_RNA.values

#cord_ATAC=pd.DataFrame(data=adata_ATAC.obs, index=adata_ATAC.obs_names, columns=['array_col','array_row'])
#adata_ATAC.obsm['spatial'] = cord_ATAC.values

#filter non-tissue barcode
barcode = pd.read_csv('position_' + sampleId+ '.txt', sep=',')
barcode = barcode.transpose()
barcode = barcode.reset_index()
barcode = barcode.rename(columns={'index': 'array'})
position = pd.read_csv("/public/home/fqjiang_gdl/SPARC/barcode_file/SPARC_barcode_filter.csv")
filte_bc = pd.merge(barcode, position, on='array', how='inner')

adata_RNA = adata_RNA[adata_RNA.obs_names.isin(filte_bc.barcode), :]
adata_ATAC = adata_ATAC[adata_ATAC.obs_names.isin(filte_bc.barcode), :]

adata_RNA.var_names_make_unique()
adata_ATAC.var_names_make_unique()

from SpatialGlue.preprocess import preprocessing
data = preprocessing(adata_RNA, adata_ATAC, datatype='Spatial-epigenome-transcriptome')

# define model
model = SpatialGlue.SpatialGlue(data, datatype='Spatial-epigenome-transcriptome', device=device)

# train model
output = model.train()

adata = adata_RNA.copy()
adata.obsm['emb_latent_RNA'] = output['emb_latent_omics1']
adata.obsm['emb_latent_ATAC'] = output['emb_latent_omics2']
adata.obsm['SpatialGlue'] = output['SpatialGlue']
adata.obsm['alpha'] = output['alpha']
adata.obsm['alpha_RNA'] = output['alpha_omics1']
adata.obsm['alpha_ATAC'] = output['alpha_omics2']

#plot UMAP and Spatial Cluster for Combined, ATAC and RNA
from SpatialGlue.utils import clustering
for i in ['emb_latent_RNA','emb_latent_ATAC','SpatialGlue']:
  clustering(adata, key=i, add_key=i, n_clusters=8, method='mclust', use_pca=True)
  sc.pp.neighbors(adata, use_rep=i, n_neighbors=30)
  sc.tl.umap(adata) # X_umap name cannot change.
  fig, ax_list = plt.subplots(1, 2, figsize=(12, 5))
  sc.pl.umap(adata, color=i, title=i, s=60, ax=ax_list[0], show=False, save = False)
  #rename X_umap to save
  adata.obsm[i+'_umap'] = adata.obsm.pop('X_umap')
#  sc.pl.embedding(adata, basis='spatial', color=i, title=i, s=90, show=False, save = sampleId+"_SpatialGlue_Spatial_"+i+".pdf")
  sq.pl.spatial_scatter(adata, img_res_key="lowres", color=i, size=15, ax=ax_list[1], save = False)
  plt.savefig(sampleId+"_SpatialGlue_"+i+".pdf")
  plt.close()
  
#save SpatialGlue.h5ad
adata.obs.rename(columns={'_index': 'index'}, inplace=True)
adata.write(sampleId + '-SpatialGlue.h5ad')

# plotting modality weight values.
plt.rcParams['figure.figsize'] = (7,3)
df = pd.DataFrame(columns=['RNA', 'ATAC', 'label'])
df['RNA'], df['ATAC'] = adata.obsm['alpha'][:, 0], adata.obsm['alpha'][:, 1]
df['label'] = adata.obs['SpatialGlue'].values # 'mclust','leiden','louvain'
df = df.set_index('label').stack().reset_index()
df.columns = ['label_SpatialGlue', 'Modality', 'Weight value']
ax = sns.violinplot(data=df, x='label_SpatialGlue', y='Weight value', hue="Modality",
                split=True, inner="quart", linewidth=1, show=False)
ax.set_title('RNA vs ATAC')
ax.set_xlabel('SpatialGlue label')
ax.legend(bbox_to_anchor=(1.2, 1.01), loc='upper right')

plt.tight_layout(w_pad=0.05)
plt.savefig(sampleId + "-modality-weight.pdf")
plt.close()

sq.pl.spatial_scatter(adata, img_res_key="lowres", color=['nCount_RNA', 'nFeature_RNA','FRiP', 'log10_unique_fragments', 'atac_percent_mito', 
                      'percent_TSS_fragments', 'percent.mt', 'percent.rp'], size=15, save = sampleId + "_Spatial_Feature_All.pdf")

#plot spatial RNA expression raw                     
sq.pl.spatial_scatter(adata, img_res_key="lowres", color=["Cldn18","Cldn5","Ftl1","Cdh5","Col1a1","Epcam","Pecam1","Lyve1","Vwf","Akap5","Cd34",
                      "Icam1",'Casz1'], size=15, save = sampleId + "_Spatial_Marker_All.pdf")


