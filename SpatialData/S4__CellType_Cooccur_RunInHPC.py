import os
import sys
import glob
import pandas as pd
import numpy as np
import anndata as ad
import spatialdata as sd
from spatialdata_io import xenium

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import squidpy as sq


saveFolder="/projects/b1146/BharatLab/Yuanqing/Projects/Spatial/IPF"
FigOut="/projects/b1042/Yuanqing/Sptial/data/output"
folderPath="/projects/b1146/BharatLab/20240530__190352__20240530_X1_Bharat"
sampleFile=glob.glob(f'{folderPath}/*__*')
disease=["IPF4","IPF1", "IPF3L", "IPF2", "Donor1", "IPF5", "IPF3R", "Donor2"]
zarr=[saveFolder+f'/Zarr/{disease[i]}.zarr' for i in range(len(disease))]
pdata=pd.DataFrame({'sam':sampleFile, 'disease':disease, 'Zarr':zarr})
pdata.disease

disease=['Donor1','Donor2',"IPF4","IPF1", "IPF3L", "IPF2","IPF5", "IPF3R"]
readPath='/projects/b1146/BharatLab/Yuanqing/Projects/Spatial/IPF/Zarr_withCT'
allZar=[sd.read_zarr(f'{readPath}/{ifile}.zarr') for ifile in disease]
allZar

allAnn=[sda.tables['table'] for sda in allZar]
allAnn[0].obs.head(5)

addvalue=500
allAnn2=allAnn.copy()
for i in range(len(disease)):
    allAnn2[i].obs['Disease']=disease[i]
    if i>0:
        maxvalue=np.max(allAnn2[i-1].obsm['spatial'][:,0])
        print(maxvalue)
        allAnn2[i].obsm['spatial'][:,0]=allAnn2[i].obsm['spatial'][:,0]+maxvalue+addvalue
        
mergDat=ad.concat(allAnn2,join='outer',axis=0)
mergDat=mergDat[mergDat.obs.CellType.isin(['unknown'])==False]
np.max(mergDat.obsm['spatial'], axis=0)

mergDat.write(f"{saveFolder}/anndata_IPF_mergForCoo.h5ad")

sq.gr.spatial_neighbors(mergDat, coord_type="generic", delaunay=True)

###Compute centrality scores
sq.gr.centrality_scores(mergDat, cluster_key="CellType")

sq.gr.co_occurrence(mergDat, cluster_key="CellType",)

sq.gr.nhood_enrichment(mergDat, cluster_key="CellType")

mergDat.write(f"{saveFolder}/anndata_IPF_Coo_CT_Final.h5ad")


from matplotlib import rcParams
# Use a vector-friendly PDF backend
rcParams["pdf.fonttype"] = 42  # Ensures text is stored as text, not paths
rcParams["ps.fonttype"] = 42

sq.pl.nhood_enrichment(mergDat, 
                       figsize=[8,6],
                       cluster_key="CellType",
                       vmin=0,  
                       vmax=20.0,
                       show=False)

# Save the current figure as a PDF
pdf_path = f"{FigOut}/allDat_IPF_enrichment_plot.pdf"
plt.savefig(pdf_path, format="pdf", bbox_inches="tight")
