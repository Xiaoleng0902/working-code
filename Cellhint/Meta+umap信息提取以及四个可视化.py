import os
import scanpy as sc
import pandas as pd
for i  in os.listdir("/data/zhouweiwei/ST/save/cancer/"):
    # print(i)
    adata=sc.read_h5ad(os.path.join("/data/zhouweiwei/ST/save/cancer/",i,"adata_after_cellhint_cluster_marker.h5ad"))
    umap=pd.DataFrame(adata.obsm["X_umap"])
    umap.set_index(adata.obs.index, inplace=True)
    umap.to_csv(os.path.join("/data/zhouweiwei/ST/save/cancer/",i,"umap.txt"),sep="\t")
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    marker=pd.DataFrame(
        {
            group + "_" + key[:10]: result[key][group]
            for group in groups
            for key in ["names","pvals","pvals_adj","scores"]
        }
    )
    marker.to_csv(os.path.join("/data/zhouweiwei/ST/save/cancer/",i,"leiden_clusters_markers.txt"),sep="\t")
    meta=pd.DataFrame(adata.obs)
    meta.to_csv(os.path.join("/data/zhouweiwei/ST/save/cancer/",i,"meta.txt"),sep="\t")
    sc.pl.umap(adata,color="copykat_anno")
    adata.uns['copykat_anno_colors']=["#5477AF","#cd3333","#669933","grey"]
    sc.pl.umap(adata,color="copykat_anno",frameon=False,save=i+"_"+"copycat.png",title="Copykat")
    sc.pl.umap(adata,color="leiden",frameon=False,save=i+"_"+"clusters.png",title="Cluster")
    sc.pl.umap(adata,color="batch",frameon=False,save=i+"_"+"after_cellhint_batch.png",title="Batch")
for i  in os.listdir("/data/zhouweiwei/ST/save/cancer/"):
    print(i)
    adata=sc.read_h5ad(os.path.join("/data/zhouweiwei/ST/save/cancer/",i,"adata_before_cellhint.h5ad"))
    sc.pl.umap(adata,color="batch",frameon=False,save=i+"_"+"before_cellhint_batch.png",title="Batch")