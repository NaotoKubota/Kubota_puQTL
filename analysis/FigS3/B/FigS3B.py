import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
sns.set_style('ticks')

def PEERcheck(QTL):
    l = []

    for i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 35, 40, 45, 50]:
        
    
        datapath = "/home/naoto/TSSchoiceQTL/data/" + QTL + "/txt/nominal_PEER"+ str(i) + "_000001_full.txt"
        
        nominal_df = pd.read_csv(datapath,
                                sep = " ", names = ["phe_id", "phe_chr", "phe_from", "phe_to", "phe_strd",
                                                   "n_var_in_cis", "dist_phe_var", "var_id", "var_chr", "var_from", "var_to",
                                                    "nom_pval", "r_squared", "slope", "best_hit"])

        # Number of best hit QTL
        l.append(nominal_df[nominal_df['best_hit'] == 1].shape[0])
        
    QTLnum_df = pd.DataFrame(l, columns = ['Number'],
                         index = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 35, 40, 45, 50]).reset_index()
    QTLnum_df.columns = ['PEER', 'Number']
    print(QTLnum_df)
    
    # make plot
    fig = plt.figure(figsize = (4, 4))
    ax = fig.add_subplot(111)
    ax = plt.plot('PEER', 'Number', data = QTLnum_df, linestyle='-', marker='o', color = '#455681')
    plt.xticks(np.arange(0, 55, 5))
    plt.xlabel('# of PEER', size = 12)
    ylabel = "# of " + QTL
    plt.ylabel(ylabel, size = 12)
    sns.despine()
    plt.tight_layout()

    pngpath = "/home/naoto/TSSchoiceQTL/png_svg/PEERcheck_" + QTL + ".png"
    plt.savefig(pngpath, dpi = 800)

PEERcheck("puQTL")