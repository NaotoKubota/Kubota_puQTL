import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')
import numpy as np

puQTL_conditional_df = pd.read_csv("/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_PEER25_01_full.txt.gz", sep = " ", header = None)
puQTL_conditional_best_df = puQTL_conditional_df[puQTL_conditional_df[22] == 1]
puQTL_conditional_best_df = puQTL_conditional_best_df[[0, 7, 8, 9, 10]]
puQTL_conditional_best_df.columns = ["prmtr", "SNP", "SNP_chr", "SNP_start", "SNP_end"]

promoter_gene_df = pd.read_csv("/home/naoto/TSSchoiceQTL/data/proActiv/promoter_gene.txt", sep = "\t",
                              names = ["prmtr_chr", "prmtr_start", "prmtr_end", "prmtr", "gene"])
singlepromoter_gene_df = promoter_gene_df.drop_duplicates(subset = "gene", keep = False)
singlepromoter_gene_df['group'] = "single"
multipromoter_gene_df = promoter_gene_df[promoter_gene_df.duplicated(subset = "gene", keep = False)]
multipromoter_gene_df['group'] = "multi"
promoter_gene_df = pd.concat([singlepromoter_gene_df, multipromoter_gene_df]).sort_values("gene")

puQTL_conditional_best_gene_df = pd.merge(puQTL_conditional_best_df, promoter_gene_df, on = "prmtr")

puQTL_group_count_df = puQTL_conditional_best_gene_df[["gene", "group"]].drop_duplicates().groupby("group").count().reset_index()

# All analyzed promoter
analyzed_promoter_df = pd.read_csv("/home/naoto/TSSchoiceQTL/data/puQTL/bed/promoter.bed", sep = "\t",
                                  names = ["prmtr_chr", "prmtr_start", "prmtr_end", "prmtr", "gene", "strand"])
analyzed_promoter_df = analyzed_promoter_df[["prmtr"]]
analyzed_promoter_df = pd.merge(analyzed_promoter_df, promoter_gene_df, on = "prmtr")
analyzed_promoter_group_df = analyzed_promoter_df[["gene", "group"]].drop_duplicates().groupby("group").count().reset_index()

# Calculate distance between puQTL and promoters
puQTL_targetgenes_df = pd.merge(puQTL_conditional_best_gene_df[puQTL_conditional_best_gene_df['group'] == "multi"], multipromoter_gene_df, on = "gene")
puQTL_targetgenes_df = puQTL_targetgenes_df[["SNP", "SNP_chr", "SNP_start", "SNP_end", "gene", "prmtr_y", "prmtr_chr_y", "prmtr_start_y", "prmtr_end_y", "group_y"]].drop_duplicates()
puQTL_targetgenes_df.columns = ["SNP", "SNP_chr", "SNP_start", "SNP_end", "gene", "prmtr", "prmtr_chr", "prmtr_start", "prmtr_end", "group"]
puQTL_targetgenes_df['distance'] = puQTL_targetgenes_df['prmtr_end'] - puQTL_targetgenes_df['SNP_end']
puQTL_targetgenes_df['abs distance'] = abs(puQTL_targetgenes_df['distance'])

puQTL_distance_df = puQTL_targetgenes_df[["SNP", "prmtr", "gene", "abs distance"]]
puQTL_distance_df["SNP-prmtr"] = puQTL_distance_df["SNP"] + "-" + puQTL_distance_df["prmtr"]

# Rank of distance
puQTL_distance_df['rank'] = puQTL_distance_df.groupby(["gene"])['abs distance'].rank()

# Merge with puQTL pairs
puQTL_SNP_prmtr_df = puQTL_conditional_best_gene_df[puQTL_conditional_best_gene_df['group'] == 'multi'][["SNP", "prmtr"]].drop_duplicates()
puQTL_SNP_prmtr_df['SNP-prmtr'] = puQTL_SNP_prmtr_df["SNP"] + "-" + puQTL_SNP_prmtr_df["prmtr"]

puQTL_SNP_prmtr_distance_rank_df = pd.merge(puQTL_SNP_prmtr_df, puQTL_distance_df, on = "SNP-prmtr")
puQTL_SNP_prmtr_distance_rank_df = puQTL_SNP_prmtr_distance_rank_df[["SNP-prmtr", "gene", "abs distance", "rank"]].drop_duplicates()

for index, row in puQTL_SNP_prmtr_distance_rank_df.iterrows():
    if row['rank'] == 1:
        puQTL_SNP_prmtr_distance_rank_df.at[index, "group"] = "Closest"
    elif row['rank'] == 2:
        puQTL_SNP_prmtr_distance_rank_df.at[index, "group"] = "2nd"
    elif row['rank'] == 3:
        puQTL_SNP_prmtr_distance_rank_df.at[index, "group"] = "3rd"
    else:
        puQTL_SNP_prmtr_distance_rank_df.at[index, "group"] = ">4th"

fig = plt.figure(figsize = (4, 3))
ax = fig.add_subplot(111)

sns.barplot(data = puQTL_SNP_prmtr_distance_rank_df.groupby("group").count().reset_index(), x = "group", y ="gene",
           palette = "Set2", order = ["Closest", "2nd", "3rd", ">4th"], axes = ax)

plt.xticks(np.arange(0, 4, 1), ["Closest", "2nd", "3rd", ">4th"])
plt.xlabel(None)
plt.ylabel("Count")

sns.despine()
plt.tight_layout()

plt.savefig("/home/naoto/TSSchoiceQTL/png_svg/puQTL_distanceToPromoter_group.png", dpi = 400)
