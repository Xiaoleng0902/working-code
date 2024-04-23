import os
import re
# import  cellhint
import pandas as pd
import scanpy as sc
import cellhint


class StIntegration():
    def __init__(self, root_path, copykat_path, save_path):
        self.root_path = root_path
        self.copykat = copykat_path
        self.save_path = save_path

    def ST_count_txt_path(self, cancer):
        all_path = []
        for i in os.listdir(self.root_path):
            cancer_type = re.findall(r"\D*", i)[0]
            if cancer_type == cancer:
                slices = os.listdir(os.path.join(self.root_path, i))
                slices_path = os.path.join(self.root_path, i)
                for slice in slices:
                    if re.findall(r"\.*ST_count", slice):
                        path = os.path.join(slices_path, slice)
                        all_path.append(path)
        return all_path

    def Copykat_patg(self, cancer):
        all_path = []
        cancer_slices = []
        for i in os.listdir(self.copykat):
            slices = os.listdir(os.path.join(self.copykat, i))
            slices_path = os.path.join(self.copykat, i)
            cancer_type = re.findall(r"\D*", i)[0]
            if cancer_type == cancer:
                for slice in slices:
                    if re.findall(r"\.*_BdyTumorCore", slice):
                        path = os.path.join(slices_path, slice)
                        a = i + "_" + slice
                        name = a.split("_BdyTumorCore.txt")
                        cancer_slices.append(name[0])
                        all_path.append(path)
        return all_path, cancer_slices

    def Read_data(self, data_path_list, copykat_path_list):
        adatas = []
        if (len(data_path_list) - 1) >= 0:
            adata0 = sc.read_text(data_path_list[0], ).T
            adata0.obs = pd.read_table(copykat_path_list[0])
            adata0.obs_names_make_unique()
            adatas.append(adata0)
        if (len(data_path_list) - 1) >= 1:
            adata1 = sc.read_text(data_path_list[1]).T
            adata1.obs = pd.read_table(copykat_path_list[1])
            adata1.obs_names_make_unique()
            adatas.append(adata1)
        if (len(data_path_list) - 1) >= 2:
            adata2 = sc.read_text(data_path_list[2]).T
            adata2.obs = pd.read_table(copykat_path_list[2])
            adata2.obs_names_make_unique()
            adatas.append(adata2)
        if (len(data_path_list) - 1) >= 3:
            adata3 = sc.read_text(data_path_list[3]).T
            adata3.obs = pd.read_table(copykat_path_list[3])
            adata3.obs_names_make_unique()
            adatas.append(adata3)
        if (len(data_path_list) - 1) >= 4:
            adata4 = sc.read_text(data_path_list[4]).T
            adata4.obs = pd.read_table(copykat_path_list[4])
            adata4.obs_names_make_unique()
            adatas.append(adata4)
        if (len(data_path_list) - 1) >= 4:
            adata5 = sc.read_text(data_path_list[5]).T
            adata5.obs = pd.read_table(copykat_path_list[5])
            adata5.obs_names_make_unique()
            adatas.append(adata5)
        if (len(data_path_list) - 1) >= 6:
            adata6 = sc.read_text(data_path_list[6]).T
            adata6.obs = pd.read_table(copykat_path_list[6])
            adata6.obs_names_make_unique()
            adatas.append(adata6)
        if (len(data_path_list) - 1) >= 7:
            adata7 = sc.read_text(data_path_list[7]).T
            adata7.obs = pd.read_table(copykat_path_list[7])
            adata7.obs_names_make_unique()
            adatas.append(adata7)
        if (len(data_path_list) - 1) >= 8:
            adata8 = sc.read_text(data_path_list[8]).T
            adata8.obs = pd.read_table(copykat_path_list[8])
            adata8.obs_names_make_unique()
            adatas.append(adata8)
        if (len(data_path_list) - 1) >= 9:
            adata9 = sc.read_text(data_path_list[9]).T
            adata9.obs = pd.read_table(copykat_path_list[9])
            adata9.obs_names_make_unique()
            adatas.append(adata9)
        if (len(data_path_list) - 1) >= 10:
            adata10 = sc.read_text(data_path_list[10]).T
            adata10.obs = pd.read_table(copykat_path_list[10])
            adata10.obs_names_make_unique()
            adatas.append(adata10)
        if (len(data_path_list) - 1) >= 11:
            adata11 = sc.read_text(data_path_list[11]).T
            adata11.obs = pd.read_table(copykat_path_list[11])
            adata11.obs_names_make_unique()
            adatas.append(adata11)
        if (len(data_path_list) - 1) >= 12:
            adata12 = sc.read_text(data_path_list[12]).T
            adata12.obs = pd.read_table(copykat_path_list[12])
            adata12.obs_names_make_unique()
            adatas.append(adata12)
        if (len(data_path_list) - 1) >= 13:
            adata13 = sc.read_text(data_path_list[13]).T
            adata13.obs = pd.read_table(copykat_path_list[13])
            adata13.obs_names_make_unique()
            adatas.append(adata13)
        if (len(data_path_list) - 1) >= 14:
            adata14 = sc.read_text(data_path_list[14]).T
            adata14.obs = pd.read_table(copykat_path_list[14])
            adata14.obs_names_make_unique()
            adatas.append(adata14)
        if (len(data_path_list) - 1) >= 15:
            adata15 = sc.read_text(data_path_list[15]).T
            adata15.obs = pd.read_table(copykat_path_list[15])
            adata15.obs_names_make_unique()
            adatas.append(adata15)
        if (len(data_path_list) - 1) >= 16:
            adata16 = sc.read_text(data_path_list[16]).T
            adata16.obs = pd.read_table(copykat_path_list[16])
            adata16.obs_names_make_unique()
            adatas.append(adata16)
        if (len(data_path_list) - 1) >= 17:
            adata17 = sc.read_text(data_path_list[17]).T
            adata17.obs = pd.read_table(copykat_path_list[17])
            adata17.obs_names_make_unique()
            adatas.append(adata17)
        if (len(data_path_list) - 1) >= 18:
            adata18 = sc.read_text(data_path_list[18]).T
            adata18.obs = pd.read_table(copykat_path_list[18])
            adata18.obs_names_make_unique()
            adatas.append(adata18)
        if (len(data_path_list) - 1) >= 19:
            adata19 = sc.read_text(data_path_list[19]).T
            adata19.obs = pd.read_table(copykat_path_list[19])
            adata19.obs_names_make_unique()
            adatas.append(adata19)
        if (len(data_path_list) - 1) >= 20:
            adata20 = sc.read_text(data_path_list[20]).T
            adata20.obs = pd.read_table(copykat_path_list[20])
            adata20.obs_names_make_unique()
            adatas.append(adata20)
        if (len(data_path_list) - 1) >= 21:
            adata21 = sc.read_text(data_path_list[21]).T
            adata21.obs = pd.read_table(copykat_path_list[21])
            adata21.obs_names_make_unique()
            adatas.append(adata21)
        if (len(data_path_list) - 1) >= 22:
            adata22 = sc.read_text(data_path_list[22]).T
            adata22.obs = pd.read_table(copykat_path_list[22])
            adata22.obs_names_make_unique()
            adatas.append(adata22)
        if (len(data_path_list) - 1) >= 23:
            adata23 = sc.read_text(data_path_list[23]).T
            adata23.obs = pd.read_table(copykat_path_list[23])
            adata23.obs_names_make_unique()
            adatas.append(adata23)
        if (len(data_path_list) - 1) >= 24:
            adata24 = sc.read_text(data_path_list[24]).T
            adata24.obs = pd.read_table(copykat_path_list[24])
            adata24.obs_names_make_unique()
            adatas.append(adata24)
        if (len(data_path_list) - 1) >= 25:
            adata25 = sc.read_text(data_path_list[25]).T
            adata25.obs = pd.read_table(copykat_path_list[25])
            adata25.obs_names_make_unique()
            adatas.append(adata25)
        if (len(data_path_list) - 1) >= 26:
            adata26 = sc.read_text(data_path_list[26]).T
            adata26.obs = pd.read_table(copykat_path_list[26])
            adata26.obs_names_make_unique()
            adatas.append(adata26)
        if (len(data_path_list) - 1) >= 27:
            adata27 = sc.read_text(data_path_list[27]).T
            adata27.obs_names_make_unique()
            adata27.obs = pd.read_table(copykat_path_list[27])
            adatas.append(adata27)
        if (len(data_path_list) - 1) >= 28:
            adata28 = sc.read_text(data_path_list[28]).T
            adata28.obs_names_make_unique()
            adata28.obs = pd.read_table(copykat_path_list[28])
            adatas.append(adata28)
        if (len(data_path_list) - 1) >= 29:
            adata29 = sc.read_text(data_path_list[29]).T
            adata29.obs_names_make_unique()
            adata29.obs = pd.read_table(copykat_path_list[29])
            adatas.append(adata29)
        if (len(data_path_list) - 1) >= 30:
            adata30 = sc.read_text(data_path_list[30]).T
            adata30.obs = pd.read_table(copykat_path_list[30])
            adata30.obs_names_make_unique()
            adatas.append(adata30)
        if (len(data_path_list) - 1) >= 31:
            adata31 = sc.read_text(data_path_list[31]).T
            adata31.obs = pd.read_table(copykat_path_list[31])
            adata31.obs_names_make_unique()
            adatas.append(adata31)
        if (len(data_path_list) - 1) >= 32:
            adata32 = sc.read_text(data_path_list[32]).T
            adata32.obs = pd.read_table(copykat_path_list[32])
            adata32.obs_names_make_unique()
            adatas.append(adata32)
        if (len(data_path_list) - 1) >= 33:
            adata33 = sc.read_text(data_path_list[33]).T
            adata33.obs = pd.read_table(copykat_path_list[33])
            adata33.obs_names_make_unique()
            adatas.append(adata33)
        if (len(data_path_list) - 1) >= 34:
            adata34 = sc.read_text(data_path_list[34]).T
            adata34.obs = pd.read_table(copykat_path_list[34])
            adata34.obs_names_make_unique()
            adatas.append(adata34)
        if (len(data_path_list) - 1) >= 35:
            adata35 = sc.read_text(data_path_list[35]).T
            adata35.obs = pd.read_table(copykat_path_list[35])
            adata35.obs_names_make_unique()
            adatas.append(adata35)
        if (len(data_path_list) - 1) >= 36:
            adata36 = sc.read_text(data_path_list[36]).T
            adata36.obs = pd.read_table(copykat_path_list[36])
            adata36.obs_names_make_unique()
            adatas.append(adata36)
        if (len(data_path_list) - 1) >= 37:
            adata37 = sc.read_text(data_path_list[37]).T
            adata37.obs_names_make_unique()
            adata37.obs = pd.read_table(copykat_path_list[37])
            adatas.append(adata37)
        if (len(data_path_list) - 1) >= 38:
            adata38 = sc.read_text(data_path_list[38]).T
            adata38.obs_names_make_unique()
            adata38.obs = pd.read_table(copykat_path_list[38])
            adatas.append(adata38)
        if (len(data_path_list) - 1) >= 39:
            adata39 = sc.read_text(data_path_list[39]).T
            adata39.obs_names_make_unique()
            adata39.obs = pd.read_table(copykat_path_list[39])
            adatas.append(adata39)
        if (len(data_path_list) - 1) >= 40:
            adata40 = sc.read_text(data_path_list[40]).T
            adata40.obs = pd.read_table(copykat_path_list[40])
            adata40.obs_names_make_unique()
            adatas.append(adata40)
        if (len(data_path_list) - 1) >= 41:
            adata41 = sc.read_text(data_path_list[41]).T
            adata41.obs = pd.read_table(copykat_path_list[41])
            adata41.obs_names_make_unique()
            adatas.append(adata41)
        if (len(data_path_list) - 1) >= 42:
            adata42 = sc.read_text(data_path_list[42]).T
            adata42.obs = pd.read_table(copykat_path_list[42])
            adata42.obs_names_make_unique()
            adatas.append(adata42)
        if (len(data_path_list) - 1) >= 43:
            adata43 = sc.read_text(data_path_list[43]).T
            adata43.obs = pd.read_table(copykat_path_list[43])
            adata43.obs_names_make_unique()
            adatas.append(adata43)
        if (len(data_path_list) - 1) >= 44:
            adata44 = sc.read_text(data_path_list[44]).T
            adata44.obs = pd.read_table(copykat_path_list[44])
            adata44.obs_names_make_unique()
            adatas.append(adata44)
        if (len(data_path_list) - 1) >= 45:
            adata45 = sc.read_text(data_path_list[45]).T
            adata45.obs = pd.read_table(copykat_path_list[45])
            adata45.obs_names_make_unique()
            adatas.append(adata45)

        adata = adatas[0].concatenate(adatas[1:], fill_value=0, join="outer")

        # adata = sc.AnnData.concatenate(adata0, adata1, adata2, adata3,)
        #                                # adata10, adata11, adata12, adata13, adata14, adata15, adata16, adata17, adata18,
        # adata19,)
        #                                adata20, adata21, adata22, adata23, adata24, adata25,)
        return adata

    def St_data_integration(self, adata):
        adata.write(os.path.join(self.save_path, "adata_concatenate.h5ad"))
        adata.layers["counts"] = adata.X.copy()
        adata.obs["batch"].value_counts()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.raw = adata
        sc.pp.highly_variable_genes(adata, batch_key='batch', subset=False, n_top_genes=3000)
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        adata.write_h5ad(os.path.join(self.save_path, "adata_before_cellhint.h5ad"))
        cellhint.integrate(adata, batch='batch')
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=0.8)
        sc.settings.verbosity = 2
        sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

        return adata


if __name__ == '__main__':
    cancers = ["brca", "cesc", "crc", "cscc", "gbm", "gist", "hgsc", "hn-as", "ipnm", "lihc", "luad", "mibc", "ovca",
               "oscc", "pcnsl", "prad", "pdac", "rcc", "skcm"]
    for cancer in cancers:
        save_path = "/data/zhouweiwei/ST/save/cancer/" + cancer.upper()
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        St = StIntegration(root_path=r"/data/zhouweiwei/ST/1ROGUE", copykat_path=r"/data/zhouweiwei/ST/copykat",
                           save_path=save_path)
        copy_path_list, cancer_slice = St.Copykat_patg(cancer)
        path_list = sorted(St.ST_count_txt_path(cancer))
        print(len(path_list))
        copy_path_list = sorted(copy_path_list)
        adata = St.Read_data(path_list, copy_path_list)
        adata1 = adata
        adata = St.St_data_integration(adata)
        slices = []
        cancer_slice = sorted(cancer_slice)
        df = pd.DataFrame(adata.obs["batch"])
        # print(df.value_counts())
        for i in df["batch"].values:
            slices.append(cancer_slice[int(i)])
        adata.obs["slice"] = pd.DataFrame(slices, columns=["slice"]).values
        adata.write(os.path.join(save_path, "adata_after_cellhint_cluster_marker.h5ad"))