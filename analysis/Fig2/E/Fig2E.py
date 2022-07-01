import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')
import numpy as np

# ChromHMM
eQTL_fenrich_chromHMM_df = pd.read_csv("/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_enrichment_chromHMM.txt", sep = "\t",
                             names = ["observed", "total", "mean expected", "sd expected", "empirical p", "lower 95% CI", "odds ratio", "upper 95% CI", "marks"])

eQTL_fenrich_chromHMM_df['empirical p log'] = -np.log10(eQTL_fenrich_chromHMM_df['empirical p'])

eQTL_fenrich_chromHMM_df.loc[eQTL_fenrich_chromHMM_df['empirical p'] < 0.05, "color"] = "red"
eQTL_fenrich_chromHMM_df.loc[~(eQTL_fenrich_chromHMM_df['empirical p'] < 0.05), "color"] = "black"

fig = plt.figure(figsize=(4, 8))
ax = fig.add_subplot(111)

# 95% CI
for NUM in list(eQTL_fenrich_chromHMM_df.index):
    ax = plt.plot([eQTL_fenrich_chromHMM_df.at[NUM, "upper 95% CI"], eQTL_fenrich_chromHMM_df.at[NUM, "lower 95% CI"]], [NUM+1, NUM+1], color = "black")

# vertical line
ax = plt.axvline(x = 1.0, color = "black", linestyle = "--", )

# scatter plot
ax = plt.scatter(data = eQTL_fenrich_chromHMM_df, x = "odds ratio", y = "marks", s = 100,
                 linewidths = 0.5, color = eQTL_fenrich_chromHMM_df["color"], edgecolors = eQTL_fenrich_chromHMM_df["color"],
                zorder = 100)
ax = plt.yticks(np.arange(1, 26 ,1), 
                ["1_TssA", "2_PromU", "3_PromD1","4_PromD2", "5_Tx5'", "6_Tx", "7_Tx3'", "8_TxWk", "9_TxReg", "10_TxEnh5'", 
                 "11_TxEnh3'", "12_TxEnhW", "13_EnhA1", "14_EnhA2", "15_EnhAF", "16_EnhW1", "17_EnhW2", "18_EnhAc", "19_DNase", "20_ZNF/Rpts", "21_Het", "22_PromP", "23_PromBiv", "24_ReprPC", "25_Quies"])

sns.despine()
plt.tight_layout()

plt.savefig("/home/naoto/TSSchoiceQTL/png_svg/fenrich_chromHMM_eQTL.png", dpi = 400)