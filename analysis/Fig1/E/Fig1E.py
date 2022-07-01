import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')
from matplotlib_venn import venn2_circles
from matplotlib import pyplot

# Import
eQTL_df = pd.read_csv("/home/naoto/TSSchoiceQTL/data/eQTL/SNP-SV/txt/conditional_PEER50_01_full.txt.gz",
                            sep = " ", names = ["phe_id", "phe_chr", "phe_from", "phe_to", "phe_strd",
                                               "n_var_in_cis", "dist_phe_var", "var_id", "var_chr", "var_from", "var_to", "rank",
                                                "fwd_pval", "fwd_r_squared", "fwd_slope", "fwd_slope_se", "fwd_best_hit", "fwd_sig",
                                                "bwd_pval", "bwd_r_squared", "bwd_slope", "bwd_slope_se", "bwd_best_hit", "bwd_sig"])
puQTL_df = pd.read_csv("/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/txt/conditional_PEER25_01_full.txt.gz",
                            sep = " ", names = ["phe_id", "phe_chr", "phe_from", "phe_to", "phe_strd",
                                               "n_var_in_cis", "dist_phe_var", "var_id", "var_chr", "var_from", "var_to", "rank",
                                                "fwd_pval", "fwd_r_squared", "fwd_slope", "fwd_slope_se", "fwd_best_hit", "fwd_sig",
                                                "bwd_pval", "bwd_r_squared", "bwd_slope", "bwd_slope_se", "bwd_best_hit", "bwd_sig"])
Promoter_Gene_df = pd.read_csv("/home/naoto/TSSchoiceQTL/data/puQTL/bed/promoter.bed",
                              sep = "\t", names = ["chr", "start", "end", "prmtr", "gene", "strand"])

puQTL_df = pd.merge(puQTL_df, Promoter_Gene_df, left_on = "phe_id", right_on = "prmtr")


eGenes = set(eQTL_df["phe_id"])
puGenes = set(puQTL_df["gene"])

A = puGenes - eGenes
B = eGenes - puGenes
C = puGenes - A
print(len(A), len(B), len(C))

g = venn2_circles((len(A), len(B), len(C)), color = "#455681")
plt.tight_layout()
plt.savefig("/home/naoto/TSSchoiceQTL/png_svg/eGenes_vs_puGenes.png", dpi = 400)
