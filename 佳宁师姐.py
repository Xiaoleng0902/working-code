import  os
import numpy as np
import  pandas as pd
from  sklearn.decomposition import  NMF
filePath="/boot3/xugang/Work/ljn/nmf/data"
data_path=[os.path.join(filePath,i) for i in os.listdir(filePath)]
# datapath=["/Users/lengxiaomifeng/Desktop/GSM5628163.txt"]
for path in data_path:
    print(f"开始读{path}路径的文件")
    df=pd.read_table(path)
    print(df.head())
    k_values=[4,5,6,7,8,9]
    genes=df.columns
    all_H=pd .DataFrame()
    for k in k_values:
        model=NMF(n_components = k,init ="nndsvd")
        W =model.fit_transform(df)
        H =model.components_
        H_with_genes =pd.DataFrame(H,columns = genes,index = [f'k{k}-{i+1}' for i in range(k)])
        all_H =pd.concat([all_H,H_with_genes], axis=0)
        print(all_H)
    all_H.to_csv(path+"_H_matrix_nndsvd.txt",sep="\t")
