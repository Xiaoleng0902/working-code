import re
import umap
import scanpy as sc
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import entropy
adata=sc.read_h5ad("/data/zhouweiwei/ST/save/outer_total/outer_total.h5ad")
sc.pl.umap(adata,color="cancer",frameon=False,title="Cancer")
colors=pd.read_table("/data/zhouweiwei/ST/save/outer_total/cancer_color.txt",sep="\t",index_col=False)
colors = colors.sort_values(by="cancer")
adata.uns['cancer_colors']=colors["color"]
sc.pl.umap(adata,color="cancer",frameon=False,title="Cancer",save="no_fill_legend_cancer_umap.pdf",legend_loc=None)
sc.pl.umap(adata,color="cancer",frameon=False,title="Cancer",legend_loc=None)