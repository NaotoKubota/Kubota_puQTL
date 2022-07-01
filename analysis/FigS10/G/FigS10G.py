import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')
import numpy as np
import scipy.stats as stats

def makeplot_eQTL(SNP, TRAIT):
    
    #### Genotype ####
    
    # Load header
    genotype_header = pd.read_csv("/home/naoto/TSSchoiceQTL/data/tsv/genotype_header_438.tsv", sep = "\t", header = None)[0]
    
    # Load genotype of SNP of interest
    SNP_df = pd.read_csv("/mnt/data5/naoto/GEUVADIS/vcf/" + SNP + ".genotypes.txt", sep = "\t", names = genotype_header)
    
    # Allele
    REF = SNP_df["REF"][0]
    ALT = SNP_df["ALT"][0]
    
    # Genotype
    REFREF = REF + "/" + REF
    REFALT = REF + "/" + ALT
    ALTALT = ALT + "/" + ALT
    
    SNP_df = SNP_df.iloc[:, 9:].T
    SNP_df.columns = [SNP]
    SNP_df = SNP_df.replace("0|0", REFREF)
    SNP_df = SNP_df.replace("0|1", REFALT)
    SNP_df = SNP_df.replace("1|0", REFALT)
    SNP_df = SNP_df.replace("1|1", ALTALT)
    
    
    #### Phenotype ####
    
    # Load trait of SNP of interest
    TRAIT_df = pd.read_csv("/home/naoto/TSSchoiceQTL/data/eQTL/bed/TMM.bed.gz", sep = "\t")
    TRAIT_df = TRAIT_df[TRAIT_df["gid"] == TRAIT]
    TRAIT_df = TRAIT_df.iloc[:, 6:].T
    TRAIT_df.columns = [TRAIT]
    
    # normalize
    TRAIT_df[TRAIT] = stats.zscore(TRAIT_df[TRAIT])

    
    #### Merge ####
    
    SNP_TRAIT_df = pd.merge(SNP_df, TRAIT_df, left_index = True, right_index = True)
    
    # Sample size
    
    REFREF_Size = "(" + str(SNP_TRAIT_df.groupby(SNP).count().at[REFREF, TRAIT]) + ")"
    REFALT_Size = "(" + str(SNP_TRAIT_df.groupby(SNP).count().at[REFALT, TRAIT]) + ")"
    if ALTALT in SNP_TRAIT_df.groupby(SNP).count().index:
        ALTALT_Size = "(" + str(SNP_TRAIT_df.groupby(SNP).count().at[ALTALT, TRAIT]) + ")"
    else:
        ALTALT_Size = "(0)"
    
    
    #### Plot ####
    
    fig = plt.figure(figsize=(3, 3))
    ax = fig.add_subplot(111)

    # Violin plot
    ax = sns.violinplot(data = SNP_TRAIT_df, x = SNP, y = TRAIT, palette = "Set2", cut = 0, order = [REFREF, REFALT, ALTALT])

    plt.title(SNP + "-" + TRAIT)
    plt.xlabel(None)
    plt.ylabel("Normalized score")
    plt.xticks(np.arange(0, 3, 1), [REFREF+"\n"+REFREF_Size, REFALT+"\n"+REFALT_Size, ALTALT+"\n"+ALTALT_Size])
    sns.despine()
    plt.tight_layout()

    plt.savefig("/home/naoto/TSSchoiceQTL/png_svg/eQTL_" + SNP + "_" + TRAIT + "_violin.png", dpi = 400)


makeplot_eQTL("28764_HG02059_ins", "MSTRG.21271")