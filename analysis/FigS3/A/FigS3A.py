import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')

# Import
metadata_df = pd.read_csv("/mnt/data5/naoto/GEUVADIS/metadata/igsr-geuvadis.tsv",
                         sep = "\t")

# Import
pca_genotype_ALL_df = pd.read_csv("/home/naoto/TSSchoiceQTL/data/puQTL/pca/SV.SNP.sort.pca",
                           sep = " ", index_col = 0)


# Set index
pca_genotype_ALL_index = pca_genotype_ALL_df.index
pca_genotype_ALL_index = [s.split("_")[-1] for s in pca_genotype_ALL_index]
pca_genotype_ALL_df.index = pca_genotype_ALL_index
pca_genotype_ALL_df = pca_genotype_ALL_df.T
pca_genotype_ALL_df = pd.merge(pca_genotype_ALL_df, metadata_df, left_index = True, right_on = "Sample name")

# ALL Genotype
fig = plt.figure(figsize=(5, 3))
ax = fig.add_subplot(111)

x = "PC1"
y = "PC2"

ax = sns.scatterplot(data = pca_genotype_ALL_df, x = x, y = y, hue = "Population code", linewidth = 0,
           edgecolor = 'black', alpha = 0.7, palette = "Set2")

plt.title("Genotype")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)

plt.tight_layout()
sns.despine()

filepath = "/home/naoto/TSSchoiceQTL/png_svg/genotype_ALL_pca_" + x + "_" + y + ".png"
plt.savefig(filepath, dpi = 400)