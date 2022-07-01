import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')
import numpy as np
import scipy.stats as stats

def makeplot_puQTL_combine(SNP, TRAIT, GENE):

    #### Genotype ####

    # Load header
    genotype_header = pd.read_csv("/home/naoto/TSSchoiceQTL/data/tsv/genotype_header_438.tsv", sep = "\t", header = None)[0]

    # Load genotype of SNP of interest
    genotype_simple = lambda x: str(x).split(":")[0]
    SNP_df = pd.read_csv("/mnt/data5/naoto/GEUVADIS/vcf/" + SNP + ".genotypes.txt", sep = "\t", names = genotype_header)
    SNP_df = SNP_df.applymap(genotype_simple)

    # Allele
    REF = SNP_df["REF"][0]
    ALT = SNP_df["ALT"][0]

    # Genotype
    REFREF = "REF/REF"
    REFALT = "REF/ALT"
    ALTALT = "ALT/ALT"

    SNP_df = SNP_df.iloc[:, 9:].T
    SNP_df.columns = [SNP]
    SNP_df = SNP_df.replace("0|0", REFREF)
    SNP_df = SNP_df.replace("0|1", REFALT)
    SNP_df = SNP_df.replace("1|0", REFALT)
    SNP_df = SNP_df.replace("1|1", ALTALT)


    #### Phenotype ####

    # Load trait of SNP of interest
    TRAIT_df = pd.read_csv("/home/naoto/TSSchoiceQTL/data/puQTL/bed/promoter_activity.bed.gz", sep = "\t")

    TRAIT_df = TRAIT_df[TRAIT_df["prmtr.promoterId"].isin(TRAIT)]
    TRAIT_list = TRAIT_df['prmtr.promoterId']
    TRAIT_df = TRAIT_df.iloc[:, 6:].T
    TRAIT_df.columns = TRAIT_list
    TRAIT_df = TRAIT_df.reset_index().melt(id_vars = "index")

    # normalize
    TRAIT_df["value"] = stats.zscore(TRAIT_df["value"])

    TRAIT_df = TRAIT_df.set_index("index")


    #### Merge ####

    SNP_TRAIT_df = pd.merge(SNP_df, TRAIT_df, left_index = True, right_index = True)

    # Sample size

    REFREF_Size = "(" + str(int((SNP_TRAIT_df.groupby(SNP).count().at[REFREF, "value"]) / len(TRAIT))) + ")"
    REFALT_Size = "(" + str(int((SNP_TRAIT_df.groupby(SNP).count().at[REFALT, "value"]) / len(TRAIT))) + ")"
    if ALTALT in SNP_TRAIT_df.groupby(SNP).count().index:
        ALTALT_Size = "(" + str(int((SNP_TRAIT_df.groupby(SNP).count().at[ALTALT, "value"]) / len(TRAIT))) + ")"
    else:
        ALTALT_Size = "(0)"

    #### Plot ####

    fig = plt.figure(figsize=(6, 3))
    ax = fig.add_subplot(111)

    # Violin plot
    ax = sns.violinplot(data = SNP_TRAIT_df, x = SNP, y = "value", hue = "prmtr.promoterId",
                        palette = "Set2", cut = 0, order = [REFREF, REFALT, ALTALT])

    plt.title(SNP + "-" + GENE)
    plt.xlabel(None)
    plt.ylabel("Normalized score")
    plt.xticks(np.arange(0, 3, 1), [REFREF+"\n"+REFREF_Size, REFALT+"\n"+REFALT_Size, ALTALT+"\n"+ALTALT_Size])
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=10)
    sns.despine()
    plt.tight_layout()

    plt.savefig("/home/naoto/TSSchoiceQTL/png_svg/puQTL_" + SNP + "_" + GENE + "_violin.png", dpi = 400)


makeplot_puQTL_combine("28764_HG02059_ins", 
                       ["prmtr.70068", "prmtr.70069"], 
                       "RSPH1")