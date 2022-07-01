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



# 1. single or multiple promoters
puQTL_conditional_SingleMulti_df = \
pd.merge(
puQTL_conditional_df[
    puQTL_conditional_df[23] == 1
],
promoter_gene_df,
left_on = 0,
right_on = "prmtr")

puQTL_conditional_SingleMulti_df['SNP-gene'] = \
puQTL_conditional_SingleMulti_df[7] + "-" + puQTL_conditional_SingleMulti_df['gene']

GROUP1 = puQTL_conditional_SingleMulti_df[["gene", "group"]].drop_duplicates().groupby('group').count().iloc[1, 0]

# 2. associated with one or multiple promoters
puQTL_conditional_SingleMulti_association_df = \
puQTL_conditional_SingleMulti_df[
    puQTL_conditional_SingleMulti_df['group'] == "multi"
][["prmtr", "SNP-gene"]].drop_duplicates().groupby("SNP-gene").count().reset_index()
puQTL_conditional_SingleMulti_association_df.columns = ["SNP-gene", "promoter count"]

puQTL_conditional_SingleMulti_promotercount_df = \
pd.merge(puQTL_conditional_SingleMulti_df, 
         puQTL_conditional_SingleMulti_association_df[["SNP-gene", "promoter count"]], 
         on = "SNP-gene")

# Number of genes of which multiple promoters are associated with one puQTL
puQTL_conditional_SingleMulti_Multiassociation_df = \
puQTL_conditional_SingleMulti_association_df[
    puQTL_conditional_SingleMulti_association_df['promoter count'] != 1
]
puQTL_conditional_SingleMulti_Multiassociation_df['gene'] = \
puQTL_conditional_SingleMulti_Multiassociation_df['SNP-gene'].str.split("-", expand = True)[1]
Multiassociation_gene = list(set(puQTL_conditional_SingleMulti_Multiassociation_df['gene']))
print(len(Multiassociation_gene))

GROUP2 = len(set(puQTL_conditional_best_gene_df['gene'])) - GROUP1 - len(Multiassociation_gene)

# 3. Effect direction
puQTL_Multiassociation_slope_df = pd.DataFrame(columns = ["prmtr", "SNP-gene", "slope", "slope direction"])

for GENE in Multiassociation_gene:
    
    SNP = \
    puQTL_conditional_SingleMulti_Multiassociation_df[
        puQTL_conditional_SingleMulti_Multiassociation_df['gene'] == GENE
    ].iloc[0, 0].split("-")[0]
    
    tmp_df = \
    puQTL_conditional_SingleMulti_df[
        (puQTL_conditional_SingleMulti_df['gene'] == GENE) &
        (puQTL_conditional_SingleMulti_df[7] == SNP)
    ][["prmtr", "SNP-gene", 20]]

    tmp_df.columns = ["prmtr", "SNP-gene", "slope"]
    tmp_df.loc[tmp_df['slope'] > 0, "slope direction"] = "+"
    tmp_df.loc[tmp_df['slope'] < 0, "slope direction"] = "-"
    
    puQTL_Multiassociation_slope_df = \
    pd.concat([puQTL_Multiassociation_slope_df, tmp_df])

puQTL_Multiassociation_slopedirection_count_df = \
puQTL_Multiassociation_slope_df[[
    "SNP-gene", "slope direction"
]].drop_duplicates().groupby("SNP-gene").count().reset_index()

# Number of genes with opposite effect of one puQTL
GROUP3 = puQTL_Multiassociation_slopedirection_count_df[
    puQTL_Multiassociation_slopedirection_count_df['slope direction'] == 2
].shape[0]

# Number of genes with concordant effect of one puQTL
GROUP4 = puQTL_Multiassociation_slopedirection_count_df[
    puQTL_Multiassociation_slopedirection_count_df['slope direction'] == 1
].shape[0]


puQTLgenes_group_count_df = pd.DataFrame([{
    "Group-1": GROUP1,
    "Group-2": GROUP2,
    "Group-3": GROUP3,
    "Group-4": GROUP4,
}]).T.reset_index()

puQTLgenes_group_count_df.columns = ["Group", "count"]
puQTLgenes_group_count_df


fig = plt.figure(figsize = (4, 4))
ax = fig.add_subplot(111)

sns.barplot(data = puQTLgenes_group_count_df, x = "Group", y = "count", ax = ax, palette = "Set2")

plt.xlabel(None)

plt.tight_layout()
sns.despine()

plt.savefig("/home/naoto/TSSchoiceQTL/png_svg/puQTLgenes_group_count.png", dpi = 400)