import pandas as pd
pd.set_option("display.max_columns", 1000)
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')
import numpy as np
import scipy.stats as stats

def comparefinemap3(trait1_p_path, trait1_pip_path, trait2_p_path, LD_path, phenotype_ID, variant_ID, output_path):
    
    # Import tables
    trait1_p_df = pd.read_csv(trait1_p_path, sep = "\t")
    trait1_pip_df = pd.read_csv(trait1_pip_path, sep = "\t", names = ["phenotype", "rs_id", "pip", "phenotype-rs_id"])
    trait1_pip_df = trait1_pip_df[trait1_pip_df['phenotype'] == phenotype_ID]
    trait2_p_df = pd.read_csv(trait2_p_path, sep = "\t")
    LD_df = pd.read_csv(LD_path, sep = "\t", names = ["chr", "pos", "variant", "r2"])
    
    # merge tables
    main_df = pd.merge(trait1_p_df, trait1_pip_df, on = "rs_id")
    main_df = pd.merge(main_df, trait2_p_df, on = "rs_id")
    main_df = pd.merge(main_df, LD_df, left_on = "rs_id", right_on = "variant")

    # log p
    main_df["log pval trait1_p"] = -np.log10(main_df["pval_nominal_x"])
    main_df["log pval trait2_p"] = -np.log10(main_df["pval_nominal_y"])

    main_lead_df = main_df[main_df['rs_id'] == variant_ID]
    main_others_df = main_df[main_df['rs_id'] != variant_ID]
    
    # Scatter plot
    fig = plt.figure(figsize=(4, 6))
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(312)
    ax3 = plt.subplot(313)

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

    ax3.set_ylabel("-log10 (nominal P)")
    
    sns.despine()
    fig.tight_layout()

    plt.savefig(output_path)


comparefinemap3(
    "/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/coloc/coloc_puQTL_prmtr.94867.txt",
    "/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/susie/susie_pip.txt",
    "/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/coloc/coloc_eQTL_MSTRG.28723.txt",
    "/home/naoto/TSSchoiceQTL/data/txt/coloc/ld_rs11154537.txt",
    "prmtr.94867",
    "rs11154537",
    "/home/naoto/TSSchoiceQTL/eps/comparefinemap_TMEM200A.eps")