import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')
import numpy as np

eQTL_fenrich_df = pd.read_csv("/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_enrichment.txt", sep = "\t",
                             names = ["observed", "total", "mean expected", "sd expected", "empirical p", "lower 95% CI", "odds ratio", "upper 95% CI", "marks"])

eQTL_fenrich_df['empirical p log'] = -np.log10(eQTL_fenrich_df['empirical p'])

eQTL_fenrich_df.loc[eQTL_fenrich_df['empirical p'] < 0.05, "color"] = "red"
eQTL_fenrich_df.loc[~(eQTL_fenrich_df['empirical p'] < 0.05), "color"] = "black"

fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(111)

# 95% CI
for NUM in list(eQTL_fenrich_df.index):
    ax = plt.plot([eQTL_fenrich_df.at[NUM, "upper 95% CI"], eQTL_fenrich_df.at[NUM, "lower 95% CI"]], [NUM, NUM], color = "black")

# vertical line
ax = plt.axvline(x = 1.0, color = "black", linestyle = "--", )

# scatter plot
ax = plt.scatter(data = eQTL_fenrich_df, x = "odds ratio", y = "marks", s = 100,
                 linewidths=0.5, color = eQTL_fenrich_df["color"], edgecolors = eQTL_fenrich_df["color"],
                zorder = 100)

sns.despine()
plt.tight_layout()

plt.savefig("/home/naoto/TSSchoiceQTL/png_svg/fenrich_eQTL.png", dpi = 400)
