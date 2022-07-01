import pandas as pd
pd.set_option("display.max_columns", 1000)
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')
import numpy as np
import scipy.stats as stats

def comparefinemap4(trait1_p_path, trait1_pip_path, trait2_p_path, trait3_p_path, LD_path, phenotype_ID, variant_ID, output_path):
    
    # Import tables
    trait1_p_df = pd.read_csv(trait1_p_path, sep = "\t")
    trait1_pip_df = pd.read_csv(trait1_pip_path, sep = "\t", names = ["phenotype", "rs_id", "pip"])
    trait1_pip_df = trait1_pip_df[trait1_pip_df['phenotype'] == phenotype_ID]
    trait2_p_df = pd.read_csv(trait2_p_path, sep = "\t")
    trait3_p_df = pd.read_csv(trait3_p_path, sep = "\t")
    LD_df = pd.read_csv(LD_path, sep = "\t", names = ["chr", "pos", "variant", "r2"])
    
    # merge tables
    main_df = pd.merge(trait1_p_df, trait1_pip_df, on = "rs_id")
    main_df = pd.merge(main_df, trait2_p_df, on = "rs_id")
    main_df = pd.merge(main_df, trait3_p_df, on = "rs_id")
    main_df = pd.merge(main_df, LD_df, left_on = "rs_id", right_on = "variant")

    # log p
    main_df["log pval trait1_p"] = -np.log10(main_df["pval_nominal_x"])
    main_df["log pval trait2_p"] = -np.log10(main_df["pval_nominal_y"])
    main_df["log pval trait3_p"] = -np.log10(main_df["pval"])

    main_lead_df = main_df[main_df['rs_id'] == variant_ID]
    main_others_df = main_df[main_df['rs_id'] != variant_ID]
    
    # Scatter plot
    fig = plt.figure(figsize=(4, 8))
    ax1 = plt.subplot(411)
    ax2 = plt.subplot(412)
    ax3 = plt.subplot(413)
    ax4 = plt.subplot(414)

    # trait1_p
    compare1 = ax1.scatter(data = main_others_df, x = "pos", y = "log pval trait1_p",
                           c = main_others_df["r2"], cmap = "viridis", vmin = 0, vmax = 1,
                           linewidth = 0.5, edgecolor = "black")
    cbar_ax = fig.add_axes([0.8, 0.85, 0.025, 0.1])
    fig.colorbar(compare1, cax=cbar_ax)
    
    compare1 = ax1.scatter(data = main_lead_df, x = "pos", y = "log pval trait1_p",
                           c = "red", linewidth = 0.5, edgecolor = "black",
                          marker = "D", s = 50)

    ax1.set_xticklabels([])
    ax1.set_ylabel("-log10 (nominal P)")

    # trait_pip
    compare2 = ax2.scatter(data = main_others_df, x = "pos", y = "pip",
                           c = main_others_df["r2"], cmap = "viridis", vmin = 0, vmax = 1,
                           linewidth = 0.5, edgecolor = "black")
    compare2 = ax2.scatter(data = main_lead_df, x = "pos", y = "pip",
                          c = "red", linewidth = 0.5, edgecolor = "black",
                          marker = "D", s = 50)

    ax2.set_ylim(-0.1, 1.1)
    ax2.set_xticklabels([])
    ax2.set_ylabel("PIP")

    # trait2_p
    compare3 = ax3.scatter(data = main_others_df, x = "pos", y = "log pval trait2_p",
                           c = main_others_df["r2"], cmap = "viridis", vmin = 0, vmax = 1,
                           linewidth = 0.5, edgecolor = "black")
    compare3 = ax3.scatter(data = main_lead_df, x = "pos", y = "log pval trait2_p",
                           c = "red", linewidth = 0.5, edgecolor = "black",
                          marker = "D", s = 50)

    ax3.set_xticklabels([])
    ax3.set_ylabel("-log10 (nominal P)")

    # trait3_p
    compare4 = ax4.scatter(data = main_others_df, x = "pos", y = "log pval trait3_p",
                           c = main_others_df["r2"], cmap = "viridis", vmin = 0, vmax = 1,
                           linewidth = 0.5, edgecolor = "black")
    compare4 = ax4.scatter(data = main_lead_df, x = "pos", y = "log pval trait3_p",
                           c = "red", linewidth = 0.5, edgecolor = "black",
                          marker = "D", s = 50)

    ax4.set_ylabel("-log10 (P)")

    sns.despine()
    fig.tight_layout()

    plt.savefig(output_path)


comparefinemap4(
    "/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/coloc/coloc_puQTL_prmtr.64998.txt",
    "/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/susie/susie_pip.txt",
    "/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/coloc/coloc_eQTL_MSTRG.19786.txt",
    "/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/coloc/coloc_GWAS_GCST003043.txt",
    "/home/naoto/TSSchoiceQTL/data/txt/coloc/ld_rs2382817.txt",
    "prmtr.64998",
    "rs2382817",
    "/home/naoto/TSSchoiceQTL/eps/comparefinemap4_TMBIM1_GCST003043.eps")