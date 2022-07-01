import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')

# import tables
prmtr_gene_df = pd.read_csv("../data/prmtrID_geneID.txt", sep = "\t", names = ["prmtr", "gene"])

puQTL_df = pd.read_csv("../data/conditional_PEER25_01_besthit.txt.gz", sep = " ", header = None)
puQTL_df = puQTL_df[[0, 7]]
puQTL_df.columns = ["prmtr", "variant"]

eQTL_df = pd.read_csv("../data/conditional_PEER50_01_full_geneID_variantID.txt.gz", sep = " ", names = ["gene", "variant"])
eQTL_df['variant-gene'] = eQTL_df['variant'] + "-" + eQTL_df['gene']

nominal_10e5_eQTL_df = pd.read_csv("../data/nominal_PEER50_000001_full_geneID_variantID.txt.gz", sep = " ", names = ["gene", "variant"])
nominal_10e5_eQTL_df['variant-gene'] = nominal_10e5_eQTL_df['variant'] + "-" + nominal_10e5_eQTL_df['gene']
# merge
puQTL_gene_df = pd.merge(puQTL_df, prmtr_gene_df, on = "prmtr")
puQTL_gene_df['variant-gene'] = puQTL_gene_df['variant'] + "-" + puQTL_gene_df['gene']

# export
puQTL_gene_df[["variant-gene"]].drop_duplicates().to_csv("../data/puQTL_gene.txt", index = None, header = None)

# Number of best hit puQTL associations
puQTL_assoc = list(puQTL_gene_df['variant-gene'])

# Number of all significant eQTL associations
eQTL_assoc = list(set(eQTL_df['variant-gene']))

# Number of nominal eQTL associations (< 10e-5)
nominal_10e5_eQTL_assoc = list(set(nominal_10e5_eQTL_df['variant-gene']))

# puQTL vs. eQTL
puQTL_eQTL_df = pd.merge(puQTL_gene_df, eQTL_df, on = 'variant-gene').drop_duplicates()

# Number of best hit puQTL association present in all eQTL associations
puQTL_assoc_in_eQTL = list(set(puQTL_eQTL_df['variant-gene']))

# stacked bar plot
fig = plt.figure(figsize = (2, 4))
ax = fig.add_subplot(111)

height1 = len(puQTL_assoc_in_eQTL) # Number of puQTL assoc called as eQTL
height2 = len(puQTL_assoc) - len(puQTL_assoc_in_eQTL) # Number of puQTL assoc NOT called as eQTL

plt.bar(0, height1, color = "grey")
plt.bar(0, height2, bottom = height1, color = '#455681')
plt.xticks([0], [None])
plt.yticks(size = 10)

sns.despine()
plt.tight_layout()

plt.savefig("../figure/association_match_puQTL_eQTL_bar.png", dpi = 400)



puQTL_ld_df = pd.read_csv("../data/ld_besthit_PEER25.txt", sep = "\t", names = ["variant", "ld variant", "r2"])

puQTL_gene_ld_df = pd.merge(puQTL_gene_df, puQTL_ld_df, on = "variant")
puQTL_gene_ld_df["ld variant-gene"] = puQTL_gene_ld_df["ld variant"] + "-" + puQTL_gene_ld_df["gene"]

# puQTL vs. eQTL
puQTL_ld_eQTL_df = pd.merge(puQTL_gene_ld_df, eQTL_df, left_on = 'ld variant-gene', right_on = 'variant-gene').drop_duplicates()

# Number of best hit puQTL association present in all eQTL associations
puQTL_ld_assoc_in_eQTL = list(set(puQTL_ld_eQTL_df['variant-gene_x']))

# Propotion of puQTL ld also called eQTL

# stacked bar plot
fig = plt.figure(figsize = (2, 4))
ax = fig.add_subplot(111)

height1 = len(puQTL_ld_assoc_in_eQTL) # Number of puQTL assoc called as eQTL
height2 = len(puQTL_assoc) - len(puQTL_ld_assoc_in_eQTL) # Number of puQTL assoc NOT called as eQTL

plt.bar(0, height1, color = "grey")
plt.bar(0, height2, bottom = height1, color = '#455681')
plt.xticks([0], [None])
plt.yticks(size = 10)

sns.despine()
plt.tight_layout()

plt.savefig("../figure/association_match_puQTL_ld_eQTL_bar.png", dpi = 400)
