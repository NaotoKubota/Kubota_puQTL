import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')
import numpy as np

def nominalCheck(QTL):
    
    datapath = "/home/naoto/TSSchoiceQTL/data/" + QTL + "/txt/nominal_PEER20_chr22.txt.gz"
        
    nominal_chr22_df = pd.read_csv(datapath,
                                sep = " ", names = ["phe_id", "phe_chr", "phe_from", "phe_to", "phe_strd",
                                                   "n_var_in_cis", "dist_phe_var", "var_id", "var_chr", "var_from", "var_to",
                                                   "nom_pval", "r_squared", "slope", "slope_se", "best_hit"])
    
    nominal_chr22_quantile_df = nominal_chr22_df[["nom_pval"]].sort_values("nom_pval").reset_index()[["nom_pval"]].reset_index()
    nominal_chr22_quantile_df['index'] = nominal_chr22_quantile_df['index'] / len(nominal_chr22_quantile_df['index'])
    nominal_chr22_quantile_df = -np.log10(nominal_chr22_quantile_df)
    pd.options.mode.use_inf_as_na = True
    nominal_chr22_quantile_df = nominal_chr22_quantile_df.dropna()

    # plot
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111)

    ax = sns.scatterplot(x = nominal_chr22_quantile_df["index"],
                                        y = nominal_chr22_quantile_df["nom_pval"],  
                                        color = '#455681', edgecolor = None, size = 1)

    x = np.arange(0, 8.5)
    y = x
    ax.plot(x, y, color = 'red', linewidth = 1)
    plt.xlim(0, 7.5)
    plt.xlabel("Expected -log10($\it{P}$-value)")
    plt.ylabel("Observed -log10($\it{P}$-value)")
    plt.legend().remove()
    sns.despine()
    plt.tight_layout()

    figname = "/home/naoto/TSSchoiceQTL/eps/" + QTL + "_nominal_chr22_qq.eps"
    plt.savefig(figname, dpi=400)


nominalCheck("puQTL")

