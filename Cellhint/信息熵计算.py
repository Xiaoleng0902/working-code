import re
import umap
import scanpy as sc
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import entropy
adata=sc.read_h5ad("/data/zhouweiwei/ST/save/outer_total/outer_total.h5ad")
data=pd.DataFrame(adata.obs["slice"])
cancer=[]
for i in data["slice"]:
    ca=re.findall(r'\D+', i.split("-")[0])[0]
    cancer.append(ca)
adata.obs["cancer"]=cancer
data=pd.DataFrame(adata.obs["slice"])
cancer=[]
for i in data["slice"]:
    ca=re.findall(r'\D+', i.split("-")[0])[0]
    cancer.append(ca)
adata.obs["cancer"]=cancer
# 计算每个簇的样本来源的熵值
cluster_entropy_source = []
for cluster_id in adata.obs['seurat_clusters'].unique():
    cluster_data =adata.obs[adata.obs['seurat_clusters'] == cluster_id]
    source_counts = cluster_data['slice'].value_counts()
    cluster_entropy_source.append(entropy(source_counts))
# 计算每个簇的癌症来源的熵值
cluster_entropy_cancer = []
for cluster_id in adata.obs['seurat_clusters'].unique():
    cluster_data = adata.obs[adata.obs['seurat_clusters'] == cluster_id]
    cancer_counts = cluster_data['cancer'].value_counts()
    cluster_entropy_cancer.append(entropy(cancer_counts))

adata.obs["entropy_source"] = adata.obs['seurat_clusters']
adata.obs["entropy_cancer"] = adata.obs['seurat_clusters']
for i, cluster_id in enumerate(adata.obs['seurat_clusters'].unique()):
    print(
        f"Cluster {cluster_id}: Sample Source Entropy = {cluster_entropy_source[i]}, Cancer Source Entropy = {cluster_entropy_cancer[i]}")
    for num, j in enumerate(adata.obs['seurat_clusters']):
        if j == cluster_id:
            adata.obs["entropy_source"][num] = cluster_entropy_source[i]
            adata.obs["entropy_cancer"][num] = cluster_entropy_cancer[i]
#归一化
from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()
nor_entropy_cancer = scaler.fit_transform(adata.obs["entropy_cancer"].values.reshape(-1, 1))
nor_entropy_source= scaler.fit_transform(adata.obs["entropy_source"].values.reshape(-1, 1))
adata.obs["nor_entropy_source"]=nor_entropy_source.round(3)
adata.obs["nor_entropy_cancer"]=nor_entropy_cancer.round(3)
