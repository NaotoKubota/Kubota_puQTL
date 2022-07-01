import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')
import numpy as np

def distance_histgram(QTL):
    
    
    datapath = "/home/naoto/TSSchoiceQTL/data/" + QTL + "/txt/conditional_PEER20_01_full.txt.gz"
    

    conditional_df = pd.read_csv(datapath,
                            sep = " ", names = ["phe_id", "phe_chr", "phe_from", "phe_to", "phe_strd",
                                               "n_var_in_cis", "dist_phe_var", "var_id", "var_chr", "var_from", "var_to", "rank",
                                                "fwd_pval", "fwd_r_squared", "fwd_slope", "fwd_slope_se", "fwd_best_hit", "fwd_sig",
                                                "bwd_pval", "bwd_r_squared", "bwd_slope", "bwd_slope_se", "bwd_best_hit", "bwd_sig"])
    conditional_df['bwd_pval_log'] = -np.log10(conditional_df['bwd_pval'])
    
    # Best hit variant
    conditional_best_df = conditional_df[
        conditional_df['bwd_best_hit'] == 1
    ]
    
    # Draw histogram
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111)

    sns.histplot(data = conditional_best_df, x = 'dist_phe_var', ax = ax, bins = 40, color = '#455681')

    ax.set_xlim(-1000000, 1000000)
    ax.set_xlabel("Distance from TSS associated with " + QTL + " (Mb)")
    plt.xticks(np.arange(-1000000, 1500000, 500000), ["-1.0", "-0.5", "0", "0.5", "1.0"])

    sns.despine()
    plt.tight_layout()
    
    figname = "/home/naoto/TSSchoiceQTL/eps/" + QTL + "_conditional_best_distance_histogram.eps"
    plt.savefig(figname, dpi=400)


distance_histgram("eQTL")