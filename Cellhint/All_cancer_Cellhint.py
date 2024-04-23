import os
import re
import cellhint
import pandas as pd
import scanpy as sc


class StIntegration():
    def __init__(self, root_path, copykat_path):
        self.root_path = root_path
        self.copykat = copykat_path

    def ST_count_txt_path(self, cancer):
        all_path = []
        cancer_slices = []
        for i in os.listdir(self.root_path):
            cancer_type = re.findall(r"\D*", i)[0]
            if cancer_type in ["brca", "cesc", "crc", "cscc", "gbm", "gist", "hgsc", "hn-as", "ipnm", "lihc", "luad",
                               "mibc", "oscc", "ovca", "pcnsl", "pdac", "prad", 'rcc', "skcm"]:
                slices = os.listdir(os.path.join(self.root_path, i))
                slices_path = os.path.join(self.root_path, i)
                for slice in slices:
                    if re.findall(r"\.*ST_count", slice):
                        path = os.path.join(slices_path, slice)
                        a = i + "_" + slice
                        name = a.split("_ST_count.txt")
                        cancer_slices.append(name[0])
                        all_path.append(path)
        return all_path, cancer_slices

    def Copykat_patg(self, cancer):
        all_path = []
        cancer_slices = []
        for i in os.listdir(self.copykat):
            cancer_type = re.findall(r"\D*", i)[0]
            if cancer_type == cancer:
                slices = os.listdir(os.path.join(self.copykat, i))
                slices_path = os.path.join(self.copykat, i)

                for slice in slices:
                    if re.findall(r"\.*_BdyTumorCore", slice):
                        path = os.path.join(slices_path, slice)
                        a = i + "_" + slice
                        name = a.split("_BdyTumorCore.txt")
                        cancer_slices.append(name[0])
                        all_path.append(path)
        return all_path, cancer_slices

    def Read_data(self, data_path_list, copykat_path_list):
        adata0 = sc.read_text(data_path_list[0], ).T
        adata0.obs = pd.read_table(copykat_path_list[0])
        adata0.obs_names_make_unique()

        adata1 = sc.read_text(data_path_list[1]).T
        adata1.obs = pd.read_table(copykat_path_list[1])
        adata1.obs_names_make_unique()
        #
        adata2 = sc.read_text(data_path_list[2]).T
        adata2.obs = pd.read_table(copykat_path_list[2])
        adata2.obs_names_make_unique()

        adata3 = sc.read_text(data_path_list[3]).T
        adata3.obs = pd.read_table(copykat_path_list[3])
        adata3.obs_names_make_unique()
        #
        adata4 = sc.read_text(data_path_list[4]).T
        adata4.obs = pd.read_table(copykat_path_list[4])
        adata4.obs_names_make_unique()

        adata5 = sc.read_text(data_path_list[5]).T
        adata5.obs = pd.read_table(copykat_path_list[5])
        adata5.obs_names_make_unique()

        adata6 = sc.read_text(data_path_list[6]).T
        adata6.obs = pd.read_table(copykat_path_list[6])
        adata6.obs_names_make_unique()

        adata7 = sc.read_text(data_path_list[7]).T
        adata7.obs = pd.read_table(copykat_path_list[7])
        adata7.obs_names_make_unique()

        adata8 = sc.read_text(data_path_list[8]).T
        adata8.obs = pd.read_table(copykat_path_list[8])
        adata8.obs_names_make_unique()
        #
        adata9 = sc.read_text(data_path_list[9]).T
        adata9.obs = pd.read_table(copykat_path_list[9])
        adata9.obs_names_make_unique()

        adata10 = sc.read_text(data_path_list[10]).T
        adata10.obs = pd.read_table(copykat_path_list[10])
        adata10.obs_names_make_unique()

        adata11 = sc.read_text(data_path_list[11]).T
        adata11.obs = pd.read_table(copykat_path_list[11])
        adata11.obs_names_make_unique()

        adata12 = sc.read_text(data_path_list[12]).T
        adata12.obs = pd.read_table(copykat_path_list[12])
        adata12.obs_names_make_unique()

        adata13 = sc.read_text(data_path_list[13]).T
        adata13.obs = pd.read_table(copykat_path_list[13])
        adata13.obs_names_make_unique()

        adata14 = sc.read_text(data_path_list[14]).T
        adata14.obs = pd.read_table(copykat_path_list[14])
        adata14.obs_names_make_unique()

        adata15 = sc.read_text(data_path_list[15]).T
        adata15.obs = pd.read_table(copykat_path_list[15])
        adata15.obs_names_make_unique()

        adata16 = sc.read_text(data_path_list[16]).T
        adata16.obs = pd.read_table(copykat_path_list[16])
        adata16.obs_names_make_unique()

        adata17 = sc.read_text(data_path_list[17]).T
        adata17.obs = pd.read_table(copykat_path_list[17])
        adata17.obs_names_make_unique()

        adata18 = sc.read_text(data_path_list[18]).T
        adata18.obs = pd.read_table(copykat_path_list[18])
        adata18.obs_names_make_unique()

        adata19 = sc.read_text(data_path_list[19]).T
        adata19.obs = pd.read_table(copykat_path_list[19])
        adata19.obs_names_make_unique()

        adata20 = sc.read_text(data_path_list[20]).T
        adata20.obs = pd.read_table(copykat_path_list[20])
        adata20.obs_names_make_unique()

        adata21 = sc.read_text(data_path_list[21]).T
        adata21.obs = pd.read_table(copykat_path_list[21])
        adata21.obs_names_make_unique()

        adata22 = sc.read_text(data_path_list[22]).T
        adata22.obs = pd.read_table(copykat_path_list[22])
        adata22.obs_names_make_unique()

        adata23 = sc.read_text(data_path_list[23]).T
        adata23.obs = pd.read_table(copykat_path_list[23])
        adata23.obs_names_make_unique()

        adata24 = sc.read_text(data_path_list[24]).T
        adata24.obs = pd.read_table(copykat_path_list[24])
        adata24.obs_names_make_unique()

        adata25 = sc.read_text(data_path_list[25]).T
        adata25.obs = pd.read_table(copykat_path_list[25])
        adata25.obs_names_make_unique()

        adata26 = sc.read_text(data_path_list[26]).T
        adata26.obs = pd.read_table(copykat_path_list[26])
        adata26.obs_names_make_unique()

        adata27 = sc.read_text(data_path_list[27]).T
        adata27.obs_names_make_unique()
        adata27.obs = pd.read_table(copykat_path_list[27])

        adata28 = sc.read_text(data_path_list[28]).T
        adata28.obs_names_make_unique()
        adata28.obs = pd.read_table(copykat_path_list[28])

        adata29 = sc.read_text(data_path_list[29]).T
        adata29.obs_names_make_unique()
        adata29.obs = pd.read_table(copykat_path_list[29])

        adata30 = sc.read_text(data_path_list[30]).T
        adata30.obs = pd.read_table(copykat_path_list[30])
        adata30.obs_names_make_unique()

        adata31 = sc.read_text(data_path_list[31]).T
        adata31.obs = pd.read_table(copykat_path_list[31])
        adata31.obs_names_make_unique()

        adata32 = sc.read_text(data_path_list[32]).T
        adata32.obs = pd.read_table(copykat_path_list[32])
        adata32.obs_names_make_unique()

        adata33 = sc.read_text(data_path_list[33]).T
        adata33.obs = pd.read_table(copykat_path_list[33])
        adata33.obs_names_make_unique()

        adata34 = sc.read_text(data_path_list[34]).T
        adata34.obs = pd.read_table(copykat_path_list[34])
        adata34.obs_names_make_unique()

        adata35 = sc.read_text(data_path_list[35]).T
        adata35.obs = pd.read_table(copykat_path_list[35])
        adata35.obs_names_make_unique()

        adata36 = sc.read_text(data_path_list[36]).T
        adata36.obs = pd.read_table(copykat_path_list[36])
        adata36.obs_names_make_unique()

        adata37 = sc.read_text(data_path_list[37]).T
        adata37.obs_names_make_unique()
        adata37.obs = pd.read_table(copykat_path_list[37])

        adata38 = sc.read_text(data_path_list[38]).T
        adata38.obs_names_make_unique()
        adata38.obs = pd.read_table(copykat_path_list[38])

        adata39 = sc.read_text(data_path_list[39]).T
        adata39.obs_names_make_unique()
        adata39.obs = pd.read_table(copykat_path_list[39])

        adata40 = sc.read_text(data_path_list[40]).T
        adata40.obs = pd.read_table(copykat_path_list[40])
        adata40.obs_names_make_unique()

        adata41 = sc.read_text(data_path_list[41]).T
        adata41.obs = pd.read_table(copykat_path_list[41])
        adata41.obs_names_make_unique()

        adata42 = sc.read_text(data_path_list[42]).T
        adata42.obs = pd.read_table(copykat_path_list[42])
        adata42.obs_names_make_unique()

        adata43 = sc.read_text(data_path_list[43]).T
        adata43.obs = pd.read_table(copykat_path_list[43])
        adata43.obs_names_make_unique()

        adata44 = sc.read_text(data_path_list[44]).T
        adata44.obs = pd.read_table(copykat_path_list[44])
        adata44.obs_names_make_unique()

        adata45 = sc.read_text(data_path_list[45]).T
        adata45.obs = pd.read_table(copykat_path_list[45])
        adata45.obs_names_make_unique()
        #
        adata46 = sc.read_text(data_path_list[46]).T
        adata46.obs = pd.read_table(copykat_path_list[46])
        adata46.obs_names_make_unique()

        adata47 = sc.read_text(data_path_list[47]).T
        adata47.obs_names_make_unique()
        adata47.obs = pd.read_table(copykat_path_list[47])

        adata48 = sc.read_text(data_path_list[48]).T
        adata48.obs_names_make_unique()
        adata48.obs = pd.read_table(copykat_path_list[48])

        adata49 = sc.read_text(data_path_list[49]).T
        adata49.obs_names_make_unique()
        adata49.obs = pd.read_table(copykat_path_list[49])

        adata50 = sc.read_text(data_path_list[50]).T
        adata50.obs = pd.read_table(copykat_path_list[50])
        adata50.obs_names_make_unique()

        adata51 = sc.read_text(data_path_list[51]).T
        adata51.obs = pd.read_table(copykat_path_list[51])
        adata51.obs_names_make_unique()

        adata52 = sc.read_text(data_path_list[52]).T
        adata52.obs = pd.read_table(copykat_path_list[52])
        adata52.obs_names_make_unique()

        adata53 = sc.read_text(data_path_list[53]).T
        adata53.obs = pd.read_table(copykat_path_list[53])
        adata53.obs_names_make_unique()

        adata54 = sc.read_text(data_path_list[54]).T
        adata54.obs = pd.read_table(copykat_path_list[54])
        adata54.obs_names_make_unique()

        adata55 = sc.read_text(data_path_list[55]).T
        adata55.obs = pd.read_table(copykat_path_list[55])
        adata55.obs_names_make_unique()
        #
        adata56 = sc.read_text(data_path_list[56]).T
        adata56.obs = pd.read_table(copykat_path_list[56])
        adata56.obs_names_make_unique()

        adata57 = sc.read_text(data_path_list[57]).T
        adata57.obs_names_make_unique()
        adata57.obs = pd.read_table(copykat_path_list[57])

        adata58 = sc.read_text(data_path_list[58]).T
        adata58.obs_names_make_unique()
        adata58.obs = pd.read_table(copykat_path_list[58])

        adata59 = sc.read_text(data_path_list[59]).T
        adata59.obs_names_make_unique()
        adata59.obs = pd.read_table(copykat_path_list[59])

        adata60 = sc.read_text(data_path_list[60]).T
        adata60.obs = pd.read_table(copykat_path_list[60])
        adata60.obs_names_make_unique()

        adata61 = sc.read_text(data_path_list[61]).T
        adata61.obs = pd.read_table(copykat_path_list[61])
        adata61.obs_names_make_unique()

        adata62 = sc.read_text(data_path_list[62]).T
        adata62.obs = pd.read_table(copykat_path_list[62])
        adata62.obs_names_make_unique()

        adata63 = sc.read_text(data_path_list[63]).T
        adata63.obs = pd.read_table(copykat_path_list[63])
        adata63.obs_names_make_unique()

        adata64 = sc.read_text(data_path_list[64]).T
        adata64.obs = pd.read_table(copykat_path_list[64])
        adata64.obs_names_make_unique()

        adata65 = sc.read_text(data_path_list[65]).T
        adata65.obs = pd.read_table(copykat_path_list[65])
        adata65.obs_names_make_unique()
        #
        adata66 = sc.read_text(data_path_list[66]).T
        adata66.obs = pd.read_table(copykat_path_list[66])
        adata66.obs_names_make_unique()

        adata67 = sc.read_text(data_path_list[67]).T
        adata67.obs_names_make_unique()
        adata67.obs = pd.read_table(copykat_path_list[67])

        adata68 = sc.read_text(data_path_list[68]).T
        adata68.obs_names_make_unique()
        adata68.obs = pd.read_table(copykat_path_list[68])

        adata69 = sc.read_text(data_path_list[69]).T
        adata69.obs_names_make_unique()
        adata69.obs = pd.read_table(copykat_path_list[69])

        adata70 = sc.read_text(data_path_list[70]).T
        adata70.obs = pd.read_table(copykat_path_list[70])
        adata70.obs_names_make_unique()

        adata71 = sc.read_text(data_path_list[71]).T
        adata71.obs = pd.read_table(copykat_path_list[71])
        adata71.obs_names_make_unique()

        adata72 = sc.read_text(data_path_list[72]).T
        adata72.obs = pd.read_table(copykat_path_list[72])
        adata72.obs_names_make_unique()

        adata73 = sc.read_text(data_path_list[73]).T
        adata73.obs = pd.read_table(copykat_path_list[73])
        adata73.obs_names_make_unique()

        adata74 = sc.read_text(data_path_list[74]).T
        adata74.obs = pd.read_table(copykat_path_list[74])
        adata74.obs_names_make_unique()

        adata75 = sc.read_text(data_path_list[75]).T
        adata75.obs = pd.read_table(copykat_path_list[75])
        adata75.obs_names_make_unique()
        #
        adata76 = sc.read_text(data_path_list[76]).T
        adata76.obs = pd.read_table(copykat_path_list[76])
        adata76.obs_names_make_unique()

        adata77 = sc.read_text(data_path_list[77]).T
        adata77.obs_names_make_unique()
        adata77.obs = pd.read_table(copykat_path_list[77])

        adata78 = sc.read_text(data_path_list[78]).T
        adata78.obs_names_make_unique()
        adata78.obs = pd.read_table(copykat_path_list[78])

        adata79 = sc.read_text(data_path_list[79]).T
        adata79.obs_names_make_unique()
        adata79.obs = pd.read_table(copykat_path_list[79])

        adata80 = sc.read_text(data_path_list[80]).T
        adata80.obs = pd.read_table(copykat_path_list[80])
        adata80.obs_names_make_unique()

        adata81 = sc.read_text(data_path_list[81]).T
        adata81.obs = pd.read_table(copykat_path_list[81])
        adata81.obs_names_make_unique()

        adata82 = sc.read_text(data_path_list[82]).T
        adata82.obs = pd.read_table(copykat_path_list[82])
        adata82.obs_names_make_unique()

        adata83 = sc.read_text(data_path_list[83]).T
        adata83.obs = pd.read_table(copykat_path_list[83])
        adata83.obs_names_make_unique()

        adata84 = sc.read_text(data_path_list[84]).T
        adata84.obs = pd.read_table(copykat_path_list[84])
        adata84.obs_names_make_unique()

        adata85 = sc.read_text(data_path_list[85]).T
        adata85.obs = pd.read_table(copykat_path_list[85])
        adata85.obs_names_make_unique()
        #
        adata86 = sc.read_text(data_path_list[86]).T
        adata86.obs = pd.read_table(copykat_path_list[86])
        adata86.obs_names_make_unique()

        adata87 = sc.read_text(data_path_list[87]).T
        adata87.obs_names_make_unique()
        adata87.obs = pd.read_table(copykat_path_list[87])

        adata88 = sc.read_text(data_path_list[88]).T
        adata88.obs_names_make_unique()
        adata88.obs = pd.read_table(copykat_path_list[88])

        adata89 = sc.read_text(data_path_list[89]).T
        adata89.obs_names_make_unique()
        adata89.obs = pd.read_table(copykat_path_list[89])

        adata90 = sc.read_text(data_path_list[90]).T
        adata90.obs = pd.read_table(copykat_path_list[90])
        adata90.obs_names_make_unique()

        adata91 = sc.read_text(data_path_list[91]).T
        adata91.obs = pd.read_table(copykat_path_list[91])
        adata91.obs_names_make_unique()

        adata92 = sc.read_text(data_path_list[92]).T
        adata92.obs = pd.read_table(copykat_path_list[92])
        adata92.obs_names_make_unique()

        adata93 = sc.read_text(data_path_list[93]).T
        adata93.obs = pd.read_table(copykat_path_list[93])
        adata93.obs_names_make_unique()

        adata94 = sc.read_text(data_path_list[94]).T
        adata94.obs = pd.read_table(copykat_path_list[94])
        adata94.obs_names_make_unique()

        adata95 = sc.read_text(data_path_list[95]).T
        adata95.obs = pd.read_table(copykat_path_list[95])
        adata95.obs_names_make_unique()
        #
        adata96 = sc.read_text(data_path_list[96]).T
        adata96.obs = pd.read_table(copykat_path_list[96])
        adata96.obs_names_make_unique()

        adata97 = sc.read_text(data_path_list[97]).T
        adata97.obs_names_make_unique()
        adata97.obs = pd.read_table(copykat_path_list[97])

        adata98 = sc.read_text(data_path_list[98]).T
        adata98.obs_names_make_unique()
        adata98.obs = pd.read_table(copykat_path_list[98])

        adata99 = sc.read_text(data_path_list[99]).T
        adata99.obs_names_make_unique()
        adata99.obs = pd.read_table(copykat_path_list[99])

        adata100 = sc.read_text(data_path_list[100], ).T
        adata100.obs = pd.read_table(copykat_path_list[100])
        adata100.obs_names_make_unique()

        adata101 = sc.read_text(data_path_list[101]).T
        adata101.obs = pd.read_table(copykat_path_list[101])
        adata101.obs_names_make_unique()
        #
        adata102 = sc.read_text(data_path_list[102]).T
        adata102.obs = pd.read_table(copykat_path_list[102])
        adata102.obs_names_make_unique()

        adata103 = sc.read_text(data_path_list[103]).T
        adata103.obs = pd.read_table(copykat_path_list[103])
        adata103.obs_names_make_unique()
        #
        adata104 = sc.read_text(data_path_list[104]).T
        adata104.obs = pd.read_table(copykat_path_list[104])
        adata104.obs_names_make_unique()

        adata105 = sc.read_text(data_path_list[105]).T
        adata105.obs = pd.read_table(copykat_path_list[105])
        adata105.obs_names_make_unique()

        adata106 = sc.read_text(data_path_list[106]).T
        adata106.obs = pd.read_table(copykat_path_list[106])
        adata106.obs_names_make_unique()

        adata107 = sc.read_text(data_path_list[107]).T
        adata107.obs = pd.read_table(copykat_path_list[107])
        adata107.obs_names_make_unique()

        adata108 = sc.read_text(data_path_list[108]).T
        adata108.obs = pd.read_table(copykat_path_list[108])
        adata108.obs_names_make_unique()
        #
        adata109 = sc.read_text(data_path_list[109]).T
        adata109.obs = pd.read_table(copykat_path_list[109])
        adata109.obs_names_make_unique()

        adata110 = sc.read_text(data_path_list[110]).T
        adata110.obs = pd.read_table(copykat_path_list[110])
        adata110.obs_names_make_unique()

        adata111 = sc.read_text(data_path_list[111]).T
        adata111.obs = pd.read_table(copykat_path_list[111])
        adata111.obs_names_make_unique()

        adata112 = sc.read_text(data_path_list[112]).T
        adata112.obs = pd.read_table(copykat_path_list[112])
        adata112.obs_names_make_unique()

        adata113 = sc.read_text(data_path_list[113]).T
        adata113.obs = pd.read_table(copykat_path_list[113])
        adata113.obs_names_make_unique()

        adata114 = sc.read_text(data_path_list[114]).T
        adata114.obs = pd.read_table(copykat_path_list[114])
        adata114.obs_names_make_unique()

        adata115 = sc.read_text(data_path_list[115]).T
        adata115.obs = pd.read_table(copykat_path_list[115])
        adata115.obs_names_make_unique()

        adata116 = sc.read_text(data_path_list[116]).T
        adata116.obs = pd.read_table(copykat_path_list[116])
        adata116.obs_names_make_unique()

        adata117 = sc.read_text(data_path_list[117]).T
        adata117.obs = pd.read_table(copykat_path_list[117])
        adata117.obs_names_make_unique()

        adata118 = sc.read_text(data_path_list[118]).T
        adata118.obs = pd.read_table(copykat_path_list[118])
        adata118.obs_names_make_unique()

        adata119 = sc.read_text(data_path_list[119]).T
        adata119.obs = pd.read_table(copykat_path_list[119])
        adata119.obs_names_make_unique()

        adata120 = sc.read_text(data_path_list[120]).T
        adata120.obs = pd.read_table(copykat_path_list[120])
        adata120.obs_names_make_unique()

        adata121 = sc.read_text(data_path_list[121]).T
        adata121.obs = pd.read_table(copykat_path_list[121])
        adata121.obs_names_make_unique()

        adata122 = sc.read_text(data_path_list[122]).T
        adata122.obs = pd.read_table(copykat_path_list[122])
        adata122.obs_names_make_unique()

        adata123 = sc.read_text(data_path_list[123]).T
        adata123.obs = pd.read_table(copykat_path_list[123])
        adata123.obs_names_make_unique()

        adata124 = sc.read_text(data_path_list[124]).T
        adata124.obs = pd.read_table(copykat_path_list[124])
        adata124.obs_names_make_unique()

        adata125 = sc.read_text(data_path_list[125]).T
        adata125.obs = pd.read_table(copykat_path_list[125])
        adata125.obs_names_make_unique()

        adata126 = sc.read_text(data_path_list[126]).T
        adata126.obs = pd.read_table(copykat_path_list[126])
        adata126.obs_names_make_unique()

        adata127 = sc.read_text(data_path_list[127]).T
        adata127.obs_names_make_unique()
        adata127.obs = pd.read_table(copykat_path_list[127])

        adata128 = sc.read_text(data_path_list[128]).T
        adata128.obs_names_make_unique()
        adata128.obs = pd.read_table(copykat_path_list[128])

        adata129 = sc.read_text(data_path_list[129]).T
        adata129.obs_names_make_unique()
        adata129.obs = pd.read_table(copykat_path_list[129])

        adata130 = sc.read_text(data_path_list[130]).T
        adata130.obs = pd.read_table(copykat_path_list[130])
        adata130.obs_names_make_unique()

        adata131 = sc.read_text(data_path_list[131]).T
        adata131.obs = pd.read_table(copykat_path_list[131])
        adata131.obs_names_make_unique()

        adata132 = sc.read_text(data_path_list[132]).T
        adata132.obs = pd.read_table(copykat_path_list[132])
        adata132.obs_names_make_unique()

        adata133 = sc.read_text(data_path_list[133]).T
        adata133.obs = pd.read_table(copykat_path_list[133])
        adata133.obs_names_make_unique()

        adata134 = sc.read_text(data_path_list[134]).T
        adata134.obs = pd.read_table(copykat_path_list[134])
        adata134.obs_names_make_unique()

        adata135 = sc.read_text(data_path_list[135]).T
        adata135.obs = pd.read_table(copykat_path_list[135])
        adata135.obs_names_make_unique()

        adata136 = sc.read_text(data_path_list[136]).T
        adata136.obs = pd.read_table(copykat_path_list[136])
        adata136.obs_names_make_unique()

        adata137 = sc.read_text(data_path_list[137]).T
        adata137.obs_names_make_unique()
        adata137.obs = pd.read_table(copykat_path_list[137])

        adata138 = sc.read_text(data_path_list[138]).T
        adata138.obs_names_make_unique()
        adata138.obs = pd.read_table(copykat_path_list[138])

        adata139 = sc.read_text(data_path_list[139]).T
        adata139.obs_names_make_unique()
        adata139.obs = pd.read_table(copykat_path_list[139])

        adata140 = sc.read_text(data_path_list[140]).T
        adata140.obs = pd.read_table(copykat_path_list[140])
        adata140.obs_names_make_unique()

        adata141 = sc.read_text(data_path_list[141]).T
        adata141.obs = pd.read_table(copykat_path_list[141])
        adata141.obs_names_make_unique()

        adata142 = sc.read_text(data_path_list[142]).T
        adata142.obs = pd.read_table(copykat_path_list[142])
        adata142.obs_names_make_unique()

        adata143 = sc.read_text(data_path_list[143]).T
        adata143.obs = pd.read_table(copykat_path_list[143])
        adata143.obs_names_make_unique()

        adata144 = sc.read_text(data_path_list[144]).T
        adata144.obs = pd.read_table(copykat_path_list[144])
        adata144.obs_names_make_unique()

        adata145 = sc.read_text(data_path_list[145]).T
        adata145.obs = pd.read_table(copykat_path_list[145])
        adata145.obs_names_make_unique()
        #
        adata146 = sc.read_text(data_path_list[146]).T
        adata146.obs = pd.read_table(copykat_path_list[146])
        adata146.obs_names_make_unique()

        adata147 = sc.read_text(data_path_list[147]).T
        adata147.obs_names_make_unique()
        adata147.obs = pd.read_table(copykat_path_list[147])

        adata148 = sc.read_text(data_path_list[148]).T
        adata148.obs_names_make_unique()
        adata148.obs = pd.read_table(copykat_path_list[148])

        adata149 = sc.read_text(data_path_list[149]).T
        adata149.obs_names_make_unique()
        adata149.obs = pd.read_table(copykat_path_list[149])

        adata150 = sc.read_text(data_path_list[150]).T
        adata150.obs = pd.read_table(copykat_path_list[150])
        adata150.obs_names_make_unique()

        adata151 = sc.read_text(data_path_list[151]).T
        adata151.obs = pd.read_table(copykat_path_list[151])
        adata151.obs_names_make_unique()

        adata152 = sc.read_text(data_path_list[152]).T
        adata152.obs = pd.read_table(copykat_path_list[152])
        adata152.obs_names_make_unique()

        adata153 = sc.read_text(data_path_list[153]).T
        adata153.obs = pd.read_table(copykat_path_list[153])
        adata153.obs_names_make_unique()

        adata154 = sc.read_text(data_path_list[154]).T
        adata154.obs = pd.read_table(copykat_path_list[154])
        adata154.obs_names_make_unique()

        adata155 = sc.read_text(data_path_list[155]).T
        adata155.obs = pd.read_table(copykat_path_list[155])
        adata155.obs_names_make_unique()
        #
        adata156 = sc.read_text(data_path_list[156]).T
        adata156.obs = pd.read_table(copykat_path_list[156])
        adata156.obs_names_make_unique()

        adata157 = sc.read_text(data_path_list[157]).T
        adata157.obs_names_make_unique()
        adata157.obs = pd.read_table(copykat_path_list[157])

        adata158 = sc.read_text(data_path_list[158]).T
        adata158.obs_names_make_unique()
        adata158.obs = pd.read_table(copykat_path_list[158])

        adata159 = sc.read_text(data_path_list[159]).T
        adata159.obs_names_make_unique()
        adata159.obs = pd.read_table(copykat_path_list[159])

        adata160 = sc.read_text(data_path_list[160]).T
        adata160.obs = pd.read_table(copykat_path_list[160])
        adata160.obs_names_make_unique()

        adata161 = sc.read_text(data_path_list[161]).T
        adata161.obs = pd.read_table(copykat_path_list[161])
        adata161.obs_names_make_unique()

        adata162 = sc.read_text(data_path_list[162]).T
        adata162.obs = pd.read_table(copykat_path_list[162])
        adata162.obs_names_make_unique()

        adata163 = sc.read_text(data_path_list[163]).T
        adata163.obs = pd.read_table(copykat_path_list[163])
        adata163.obs_names_make_unique()

        adata164 = sc.read_text(data_path_list[164]).T
        adata164.obs = pd.read_table(copykat_path_list[164])
        adata164.obs_names_make_unique()

        adata165 = sc.read_text(data_path_list[165]).T
        adata165.obs = pd.read_table(copykat_path_list[165])
        adata165.obs_names_make_unique()
        #
        adata166 = sc.read_text(data_path_list[166]).T
        adata166.obs = pd.read_table(copykat_path_list[166])
        adata166.obs_names_make_unique()

        adata167 = sc.read_text(data_path_list[167]).T
        adata167.obs_names_make_unique()
        adata167.obs = pd.read_table(copykat_path_list[167])

        adata168 = sc.read_text(data_path_list[168]).T
        adata168.obs_names_make_unique()
        adata168.obs = pd.read_table(copykat_path_list[168])

        adata169 = sc.read_text(data_path_list[169]).T
        adata169.obs_names_make_unique()
        adata169.obs = pd.read_table(copykat_path_list[169])

        adata170 = sc.read_text(data_path_list[170]).T
        adata170.obs = pd.read_table(copykat_path_list[170])
        adata170.obs_names_make_unique()

        adata171 = sc.read_text(data_path_list[171]).T
        adata171.obs = pd.read_table(copykat_path_list[171])
        adata171.obs_names_make_unique()

        adata172 = sc.read_text(data_path_list[172]).T
        adata172.obs = pd.read_table(copykat_path_list[172])
        adata172.obs_names_make_unique()

        adata173 = sc.read_text(data_path_list[173]).T
        adata173.obs = pd.read_table(copykat_path_list[173])
        adata173.obs_names_make_unique()

        adata174 = sc.read_text(data_path_list[174]).T
        adata174.obs = pd.read_table(copykat_path_list[174])
        adata174.obs_names_make_unique()

        adata175 = sc.read_text(data_path_list[175]).T
        adata175.obs = pd.read_table(copykat_path_list[175])
        adata175.obs_names_make_unique()
        #
        adata176 = sc.read_text(data_path_list[176]).T
        adata176.obs = pd.read_table(copykat_path_list[176])
        adata176.obs_names_make_unique()

        adata177 = sc.read_text(data_path_list[177]).T
        adata177.obs_names_make_unique()
        adata177.obs = pd.read_table(copykat_path_list[177])

        adata178 = sc.read_text(data_path_list[178]).T
        adata178.obs_names_make_unique()
        adata178.obs = pd.read_table(copykat_path_list[178])

        adata179 = sc.read_text(data_path_list[179]).T
        adata179.obs_names_make_unique()
        adata179.obs = pd.read_table(copykat_path_list[179])

        adata180 = sc.read_text(data_path_list[180]).T
        adata180.obs = pd.read_table(copykat_path_list[180])
        adata180.obs_names_make_unique()

        adata181 = sc.read_text(data_path_list[181]).T
        adata181.obs = pd.read_table(copykat_path_list[181])
        adata181.obs_names_make_unique()

        adata182 = sc.read_text(data_path_list[182]).T
        adata182.obs = pd.read_table(copykat_path_list[182])
        adata182.obs_names_make_unique()

        adata183 = sc.read_text(data_path_list[183]).T
        adata183.obs = pd.read_table(copykat_path_list[183])
        adata183.obs_names_make_unique()

        adata184 = sc.read_text(data_path_list[184]).T
        adata184.obs = pd.read_table(copykat_path_list[184])
        adata184.obs_names_make_unique()

        adata185 = sc.read_text(data_path_list[185]).T
        adata185.obs = pd.read_table(copykat_path_list[185])
        adata185.obs_names_make_unique()
        #
        adata186 = sc.read_text(data_path_list[186]).T
        adata186.obs = pd.read_table(copykat_path_list[186])
        adata186.obs_names_make_unique()

        adata187 = sc.read_text(data_path_list[187]).T
        adata187.obs_names_make_unique()
        adata187.obs = pd.read_table(copykat_path_list[187])

        adata188 = sc.read_text(data_path_list[188]).T
        adata188.obs_names_make_unique()
        adata188.obs = pd.read_table(copykat_path_list[188])

        adata189 = sc.read_text(data_path_list[189]).T
        adata189.obs_names_make_unique()
        adata189.obs = pd.read_table(copykat_path_list[189])

        adata190 = sc.read_text(data_path_list[190]).T
        adata190.obs = pd.read_table(copykat_path_list[190])
        adata190.obs_names_make_unique()

        adata191 = sc.read_text(data_path_list[191]).T
        adata191.obs = pd.read_table(copykat_path_list[191])
        adata191.obs_names_make_unique()

        adata192 = sc.read_text(data_path_list[192]).T
        adata192.obs = pd.read_table(copykat_path_list[192])
        adata192.obs_names_make_unique()

        adata193 = sc.read_text(data_path_list[193]).T
        adata193.obs = pd.read_table(copykat_path_list[193])
        adata193.obs_names_make_unique()

        adata194 = sc.read_text(data_path_list[194]).T
        adata194.obs = pd.read_table(copykat_path_list[194])
        adata194.obs_names_make_unique()

        adata195 = sc.read_text(data_path_list[195]).T
        adata195.obs = pd.read_table(copykat_path_list[195])
        adata195.obs_names_make_unique()
        #
        adata196 = sc.read_text(data_path_list[196]).T
        adata196.obs = pd.read_table(copykat_path_list[196])
        adata196.obs_names_make_unique()

        adata197 = sc.read_text(data_path_list[197]).T
        adata197.obs_names_make_unique()
        adata197.obs = pd.read_table(copykat_path_list[197])

        adata198 = sc.read_text(data_path_list[198]).T
        adata198.obs_names_make_unique()
        adata198.obs = pd.read_table(copykat_path_list[198])

        adata199 = sc.read_text(data_path_list[199]).T
        adata199.obs_names_make_unique()
        adata199.obs = pd.read_table(copykat_path_list[199])

        adata200 = sc.read_text(data_path_list[200], ).T
        adata200.obs = pd.read_table(copykat_path_list[200])
        adata200.obs_names_make_unique()

        adata201 = sc.read_text(data_path_list[201]).T
        adata201.obs = pd.read_table(copykat_path_list[201])
        adata201.obs_names_make_unique()
        #
        adata202 = sc.read_text(data_path_list[202]).T
        adata202.obs = pd.read_table(copykat_path_list[202])
        adata202.obs_names_make_unique()

        adata203 = sc.read_text(data_path_list[203]).T
        adata203.obs = pd.read_table(copykat_path_list[203])
        adata203.obs_names_make_unique()
        #
        adata204 = sc.read_text(data_path_list[204]).T
        adata204.obs = pd.read_table(copykat_path_list[204])
        adata204.obs_names_make_unique()

        adata205 = sc.read_text(data_path_list[205]).T
        adata205.obs = pd.read_table(copykat_path_list[205])
        adata205.obs_names_make_unique()

        adata206 = sc.read_text(data_path_list[206]).T
        adata206.obs = pd.read_table(copykat_path_list[206])
        adata206.obs_names_make_unique()

        adata207 = sc.read_text(data_path_list[207]).T
        adata207.obs = pd.read_table(copykat_path_list[207])
        adata207.obs_names_make_unique()

        adata208 = sc.read_text(data_path_list[208]).T
        adata208.obs = pd.read_table(copykat_path_list[208])
        adata208.obs_names_make_unique()
        #
        adata209 = sc.read_text(data_path_list[209]).T
        adata209.obs = pd.read_table(copykat_path_list[209])
        adata209.obs_names_make_unique()

        adata210 = sc.read_text(data_path_list[210]).T
        adata210.obs = pd.read_table(copykat_path_list[210])
        adata210.obs_names_make_unique()

        adata211 = sc.read_text(data_path_list[211]).T
        adata211.obs = pd.read_table(copykat_path_list[211])
        adata211.obs_names_make_unique()

        adata212 = sc.read_text(data_path_list[212]).T
        adata212.obs = pd.read_table(copykat_path_list[212])
        adata212.obs_names_make_unique()

        adata213 = sc.read_text(data_path_list[213]).T
        adata213.obs = pd.read_table(copykat_path_list[213])
        adata213.obs_names_make_unique()

        adata214 = sc.read_text(data_path_list[214]).T
        adata214.obs = pd.read_table(copykat_path_list[214])
        adata214.obs_names_make_unique()

        adata215 = sc.read_text(data_path_list[215]).T
        adata215.obs = pd.read_table(copykat_path_list[215])
        adata215.obs_names_make_unique()

        adata216 = sc.read_text(data_path_list[216]).T
        adata216.obs = pd.read_table(copykat_path_list[216])
        adata216.obs_names_make_unique()

        adata217 = sc.read_text(data_path_list[217]).T
        adata217.obs = pd.read_table(copykat_path_list[217])
        adata217.obs_names_make_unique()

        adata218 = sc.read_text(data_path_list[218]).T
        adata218.obs = pd.read_table(copykat_path_list[218])
        adata218.obs_names_make_unique()

        adata219 = sc.read_text(data_path_list[219]).T
        adata219.obs = pd.read_table(copykat_path_list[219])
        adata219.obs_names_make_unique()

        adata220 = sc.read_text(data_path_list[220]).T
        adata220.obs = pd.read_table(copykat_path_list[220])
        adata220.obs_names_make_unique()

        adata221 = sc.read_text(data_path_list[221]).T
        adata221.obs = pd.read_table(copykat_path_list[221])
        adata221.obs_names_make_unique()

        adata222 = sc.read_text(data_path_list[222]).T
        adata222.obs = pd.read_table(copykat_path_list[222])
        adata222.obs_names_make_unique()

        adata223 = sc.read_text(data_path_list[223]).T
        adata223.obs = pd.read_table(copykat_path_list[223])
        adata223.obs_names_make_unique()

        adata224 = sc.read_text(data_path_list[224]).T
        adata224.obs = pd.read_table(copykat_path_list[224])
        adata224.obs_names_make_unique()

        adata225 = sc.read_text(data_path_list[225]).T
        adata225.obs = pd.read_table(copykat_path_list[225])
        adata225.obs_names_make_unique()

        adata226 = sc.read_text(data_path_list[226]).T
        adata226.obs = pd.read_table(copykat_path_list[226])
        adata226.obs_names_make_unique()

        adata227 = sc.read_text(data_path_list[227]).T
        adata227.obs_names_make_unique()
        adata227.obs = pd.read_table(copykat_path_list[227])

        adata228 = sc.read_text(data_path_list[228]).T
        adata228.obs_names_make_unique()
        adata228.obs = pd.read_table(copykat_path_list[228])

        adata = sc.AnnData.concatenate(adata0, adata1, adata2, adata3, adata4, adata5, adata6, adata7, adata8, adata9,
                                       adata10, adata11, adata12, adata13, adata14, adata15, adata16,
                                       adata17, adata18, adata19, adata20, adata21, adata22, adata23, adata24, adata25,
                                       adata26, adata27, adata28, adata29,
                                       adata30, adata31, adata32, adata33, adata34, adata35, adata36, adata37, adata38,
                                       adata39,
                                       adata40, adata41, adata42, adata43, adata44, adata45, adata46, adata47, adata48,
                                       adata49,
                                       adata50, adata51, adata52, adata53, adata54, adata55, adata56, adata57, adata58,
                                       adata59,
                                       adata60, adata61, adata62, adata63, adata64, adata65, adata66, adata67, adata68,
                                       adata69,
                                       adata70, adata71, adata72, adata73, adata74, adata75, adata76, adata77, adata78,
                                       adata79,
                                       adata80, adata81, adata82, adata83, adata84, adata85, adata86, adata87, adata88,
                                       adata89,
                                       adata90, adata91, adata92, adata93, adata94, adata95, adata96, adata97, adata98,
                                       adata99,
                                       adata100, adata101, adata102, adata103, adata104, adata105, adata106, adata107,
                                       adata108, adata109, adata110,
                                       adata111, adata112, adata113, adata114, adata115, adata116,
                                       adata117, adata118, adata119, adata120, adata121, adata122, adata123, adata124,
                                       adata125, adata126, adata127, adata128, adata129,
                                       adata130, adata131, adata132, adata133, adata134, adata135, adata136, adata137,
                                       adata138, adata139,
                                       adata140, adata141, adata142, adata143, adata144, adata145, adata146, adata147,
                                       adata148, adata149,
                                       adata150, adata151, adata152, adata153, adata154, adata155, adata156, adata157,
                                       adata158, adata159,
                                       adata160, adata161, adata162, adata163, adata164, adata165, adata166, adata167,
                                       adata168, adata169,
                                       adata170, adata171, adata172, adata173, adata174, adata175, adata176, adata177,
                                       adata178, adata179,
                                       adata180, adata181, adata182, adata183, adata184, adata185, adata186, adata187,
                                       adata188, adata189,
                                       adata190, adata191, adata192, adata193, adata194, adata195, adata196, adata197,
                                       adata198, adata199,
                                       adata200, adata201, adata202, adata203, adata204, adata205, adata206, adata207,
                                       adata208, adata209, adata210,
                                       adata211, adata212, adata213, adata214, adata215, adata216,
                                       adata217, adata218, adata219, adata220, adata221, adata222, adata223, adata224,
                                       adata225, adata226, adata227, adata228
                                       ,
                                       fill_value=0, join="outer"
                                       )

        # adata = sc.AnnData.concatenate(adata0, adata1, adata2, adata4,)
        #                                # adata10, adata11, adata12, adata13, adata14, adata19, adata16, adata17, adata18,
        # adata19,)
        #                                adata20, adata21, adata22, adata23, adata24, adata25,)
        return adata

    def St_data_integration(self, adata):
        adata.obs["batch"].value_counts()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.raw = adata
        sc.pp.highly_variable_genes(adata, batch_key='batch', subset=False, n_top_genes=3000)
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        sc.pl.umap(adata, color=['batch', ], wspace=0.5)
        cellhint.integrate(adata, batch='batch')
        sc.tl.umap(adata)
        sc.pl.umap(adata, color=['batch', ], wspace=0.5)
        return adata


if __name__ == '__main__':
    St = StIntegration(root_path=r"/data/zhouweiwei/ST/1ROGUE", copykat_path=None)

    path_list = [
        '/data/zhouweiwei/ST/1ROGUE/brca01/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca01/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca01/slice3_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca07/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca07/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca08/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca08/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca08/slice3_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca08/slice4_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca09/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca10/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca11/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca12/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca13/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca14/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca14/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca15/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca15/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca16/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca26/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca26/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca26/slice3_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca26/slice4_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca26/slice5_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca26/slice6_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca26/slice7_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/brca26/slice8_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/cesc01/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/cesc02/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/cesc02/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/cesc02/slice3_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/crc01/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/crc02/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/crc03/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/crc04/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/crc04/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/crc04/slice3_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/crc04/slice4_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/crc05/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/crc05/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/crc06/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/crc06/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/crc07/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/crc07/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/crc08/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/crc08/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/crc09/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/crc10/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/crc10/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/crc10/slice3_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/crc10/slice4_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/cscc01/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/cscc01/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/cscc01/slice3_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/cscc01/slice4_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/cscc02/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/cscc02/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/cscc02/slice3_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/cscc02/slice4_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm01/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm02/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm03/slice10_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm03/slice11_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm03/slice12_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm03/slice13_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm03/slice14_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm03/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm03/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm03/slice3_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm03/slice4_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm03/slice5_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm03/slice6_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm03/slice7_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm03/slice8_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm03/slice9_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm04/slice10_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm04/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm04/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm04/slice3_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm04/slice4_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm04/slice5_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm04/slice6_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm04/slice7_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm04/slice8_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm04/slice9_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm05/slice10_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm05/slice11_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm05/slice12_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm05/slice13_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm05/slice14_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm05/slice15_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm05/slice16_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm05/slice17_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm05/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm05/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm05/slice3_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm05/slice4_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm05/slice5_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm05/slice6_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm05/slice7_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm05/slice8_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm05/slice9_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm06/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gbm06/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gist01/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/gist01/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/hgsc01/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/hgsc01/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/hgsc01/slice3_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/hgsc01/slice4_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/hn-as01/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/hn-as01/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/hn-as01/slice3_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/hn-as01/slice4_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/hn-as01/slice5_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/hn-as01/slice6_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/hn-as01/slice7_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/hn-as01/slice8_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/hn-as02/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/hn-as02/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/hn-as02/slice3_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/hn-as02/slice4_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/hn-as03/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/hn-as03/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ipnm01/slice10_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ipnm01/slice11_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ipnm01/slice12_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ipnm01/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ipnm01/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ipnm01/slice3_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ipnm01/slice4_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ipnm01/slice5_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ipnm01/slice6_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ipnm01/slice7_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ipnm01/slice8_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ipnm01/slice9_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc01/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc02/slice10_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc02/slice11_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc02/slice12_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc02/slice13_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc02/slice14_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc02/slice15_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc02/slice16_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc02/slice17_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc02/slice18_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc02/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc02/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc02/slice3_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc02/slice4_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc02/slice5_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc02/slice6_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc02/slice7_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc02/slice8_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc02/slice9_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc03/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc03/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc03/slice3_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc03/slice4_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc03/slice5_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc03/slice6_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc03/slice7_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/lihc03/slice8_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/luad01/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/luad01/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/luad02/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/luad02/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/luad03/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/luad03/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/luad04/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/luad05/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/luad06/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/mibc01/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/mibc01/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/mibc01/slice3_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/mibc01/slice4_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/oscc01/slice10_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/oscc01/slice11_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/oscc01/slice12_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/oscc01/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/oscc01/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/oscc01/slice3_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/oscc01/slice4_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/oscc01/slice5_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/oscc01/slice6_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/oscc01/slice7_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/oscc01/slice8_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/oscc01/slice9_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ovca01/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ovca01/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ovca02/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ovca03/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ovca04/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ovca05/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ovca06/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ovca06/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ovca06/slice3_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ovca06/slice4_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ovca06/slice5_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ovca06/slice6_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ovca07/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ovca07/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ovca07/slice3_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ovca07/slice4_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ovca07/slice5_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ovca07/slice6_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ovca07/slice7_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/ovca07/slice8_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/pcnsl01/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/pcnsl01/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/pcnsl01/slice3_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/pcnsl01/slice4_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/pdac02/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/pdac03/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/pdac03/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/pdac03/slice3_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/prad02/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/prad07/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/prad08/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/prad09/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/rcc01/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/rcc01/slice2_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/rcc01/slice3_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/rcc01/slice4_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/rcc01/slice5_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/rcc02/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/rcc03/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/skcm12/slice1_ST_count.txt',
        '/data/zhouweiwei/ST/1ROGUE/skcm13/slice1_ST_count.txt']
    copy_path_list = ['/data/zhouweiwei/ST/copykat/brca01/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca01/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca01/slice3_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca07/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca07/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca08/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca08/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca08/slice3_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca08/slice4_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca09/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca10/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca11/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca12/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca13/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca14/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca14/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca15/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca15/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca16/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca26/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca26/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca26/slice3_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca26/slice4_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca26/slice5_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca26/slice6_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca26/slice7_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/brca26/slice8_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/cesc01/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/cesc02/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/cesc02/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/cesc02/slice3_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/crc01/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/crc02/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/crc03/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/crc04/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/crc04/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/crc04/slice3_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/crc04/slice4_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/crc05/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/crc05/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/crc06/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/crc06/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/crc07/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/crc07/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/crc08/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/crc08/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/crc09/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/crc10/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/crc10/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/crc10/slice3_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/crc10/slice4_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/cscc01/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/cscc01/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/cscc01/slice3_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/cscc01/slice4_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/cscc02/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/cscc02/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/cscc02/slice3_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/cscc02/slice4_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm01/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm02/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm03/slice10_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm03/slice11_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm03/slice12_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm03/slice13_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm03/slice14_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm03/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm03/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm03/slice3_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm03/slice4_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm03/slice5_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm03/slice6_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm03/slice7_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm03/slice8_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm03/slice9_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm04/slice10_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm04/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm04/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm04/slice3_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm04/slice4_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm04/slice5_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm04/slice6_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm04/slice7_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm04/slice8_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm04/slice9_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm05/slice10_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm05/slice11_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm05/slice12_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm05/slice13_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm05/slice14_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm05/slice15_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm05/slice16_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm05/slice17_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm05/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm05/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm05/slice3_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm05/slice4_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm05/slice5_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm05/slice6_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm05/slice7_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm05/slice8_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm05/slice9_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm06/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gbm06/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gist01/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/gist01/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/hgsc01/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/hgsc01/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/hgsc01/slice3_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/hgsc01/slice4_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/hn-as01/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/hn-as01/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/hn-as01/slice3_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/hn-as01/slice4_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/hn-as01/slice5_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/hn-as01/slice6_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/hn-as01/slice7_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/hn-as01/slice8_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/hn-as02/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/hn-as02/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/hn-as02/slice3_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/hn-as02/slice4_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/hn-as03/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/hn-as03/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ipnm01/slice10_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ipnm01/slice11_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ipnm01/slice12_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ipnm01/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ipnm01/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ipnm01/slice3_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ipnm01/slice4_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ipnm01/slice5_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ipnm01/slice6_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ipnm01/slice7_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ipnm01/slice8_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ipnm01/slice9_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc01/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc02/slice10_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc02/slice11_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc02/slice12_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc02/slice13_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc02/slice14_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc02/slice15_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc02/slice16_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc02/slice17_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc02/slice18_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc02/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc02/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc02/slice3_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc02/slice4_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc02/slice5_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc02/slice6_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc02/slice7_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc02/slice8_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc02/slice9_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc03/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc03/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc03/slice3_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc03/slice4_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc03/slice5_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc03/slice6_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc03/slice7_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/lihc03/slice8_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/luad01/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/luad01/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/luad02/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/luad02/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/luad03/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/luad03/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/luad04/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/luad05/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/luad06/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/mibc01/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/mibc01/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/mibc01/slice3_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/mibc01/slice4_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/oscc01/slice10_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/oscc01/slice11_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/oscc01/slice12_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/oscc01/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/oscc01/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/oscc01/slice3_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/oscc01/slice4_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/oscc01/slice5_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/oscc01/slice6_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/oscc01/slice7_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/oscc01/slice8_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/oscc01/slice9_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ovca01/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ovca01/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ovca02/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ovca03/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ovca04/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ovca05/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ovca06/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ovca06/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ovca06/slice3_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ovca06/slice4_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ovca06/slice5_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ovca06/slice6_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ovca07/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ovca07/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ovca07/slice3_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ovca07/slice4_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ovca07/slice5_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ovca07/slice6_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ovca07/slice7_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/ovca07/slice8_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/pcnsl01/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/pcnsl01/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/pcnsl01/slice3_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/pcnsl01/slice4_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/pdac02/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/pdac03/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/pdac03/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/pdac03/slice3_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/prad02/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/prad07/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/prad08/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/prad09/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/rcc01/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/rcc01/slice2_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/rcc01/slice3_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/rcc01/slice4_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/rcc01/slice5_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/rcc02/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/rcc03/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/skcm12/slice1_BdyTumorCore.txt',
                      '/data/zhouweiwei/ST/copykat/skcm13/slice1_BdyTumorCore.txt']
    adata = St.Read_data(data_path_list=path_list, copykat_path_list=copy_path_list)
    adata1 = adata
    adata = St.St_data_integration(adata)
    slices = []
    df = pd.DataFrame(adata.obs["batch"])
    # print(df.value_counts())
    all_path, cancer_slices=St.ST_count_txt_path()
    for i in df["batch"].values:
        slices.append(cancer_slices[int(i)])
    adata.obs["slice"] = pd.DataFrame(slices, columns=["slice"]).values
    adata.write(r"/data/zhouweiwei/ST/save/total.h5ad")
