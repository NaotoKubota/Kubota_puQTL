import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')
import numpy as np

puQTL_PEER25_fenrich_df = pd.read_csv("/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_PEER25_01_enrichment.txt", sep = "\t",
                             names = ["observed", "total", "mean expected", "sd expected", "empirical p", "lower 95% CI", "odds ratio", "upper 95% CI", "marks"])

puQTL_PEER25_fenrich_df['empirical p log'] = -np.log10(puQTL_PEER25_fenrich_df['empirical p'])

puQTL_PEER25_fenrich_df.loc[puQTL_PEER25_fenrich_df['empirical p'] < 0.05, "color"] = "red"
puQTL_PEER25_fenrich_df.loc[~(puQTL_PEER25_fenrich_df['empirical p'] < 0.05), "color"] = "black"

fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(111)

# 95% CI
for NUM in list(puQTL_PEER25_fenrich_df.index):
    ax = plt.plot([puQTL_PEER25_fenrich_df.at[NUM, "upper 95% CI"], puQTL_PEER25_fenrich_df.at[NUM, "lower 95% CI"]], [NUM, NUM], color = "black")

# vertical line
ax = plt.axvline(x = 1.0, color = "black", linestyle = "--", )

# scatter plot
ax = plt.scatter(data = puQTL_PEER25_fenrich_df, x = "odds ratio", y = "marks", s = 100,
                 linewidths=0.5, color = puQTL_PEER25_fenrich_df["color"], edgecolors = puQTL_PEER25_fenrich_df["color"],
                zorder = 100)

sns.despine()
plt.tight_layout()

plt.savefig("/home/naoto/TSSchoiceQTL/png_svg/fenrich_puQTL_PEER25.png", dpi = 400)