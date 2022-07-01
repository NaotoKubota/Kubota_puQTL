import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')
import numpy as np
import scipy.stats as stats


header_CAGE = ['chr', 'start', 'end', "CAGE peak", 'CAGE score', 'strand', '1st exon ID']
CAGE_GM12878_df = pd.read_csv("/home/naoto/data6/ENCODE/GM12878/CAGE/Expression_TSS_cluster_1stExon.bed", sep = "\t",
                             names = header_CAGE)
header_proActiv = ['chr', 'start', 'end', 'proActiv ID', 'ENSG', 'strand', 'internalPromoter', 'promoterPosition', 'GM12878SJ.out', 'K562SJ.out', '1st exon ID']
proActiv_GM12878_df = pd.read_csv("/home/naoto/TSSchoiceQTL/data/proActiv_merged_sorted_1stExon.bed", sep = "\t",
                             names = header_proActiv)

# Max CAGE score of 1st exon
CAGE_GM12878_1stExon_df = \
CAGE_GM12878_df.groupby('1st exon ID')[["CAGE score"]].max().reset_index()

# proActiv score of 1st exon
proActiv_GM12878_1stExon_df = \
proActiv_GM12878_df[["1st exon ID", "GM12878SJ.out", "proActiv ID", "internalPromoter"]]

# Merge
CAGEvsproActiv_GM12878_1stExon_df = pd.merge(CAGE_GM12878_1stExon_df, proActiv_GM12878_1stExon_df, on = '1st exon ID')
CAGEvsproActiv_GM12878_1stExon_df.columns = ['1st exon ID', 'CAGE', 'proActiv', "proActiv ID", "internalPromoter"]
CAGEvsproActiv_GM12878_1stExon_df.head()

# Remove inactive promoter
CAGEvsproActiv_GM12878_1stExon_active_df = \
CAGEvsproActiv_GM12878_1stExon_df[
    ~((CAGEvsproActiv_GM12878_1stExon_df['CAGE'] == 0) & \
    (CAGEvsproActiv_GM12878_1stExon_df['proActiv'] == 0))
]

CAGEvsproActiv_GM12878_1stExon_active_plot_df = \
np.log10(CAGEvsproActiv_GM12878_1stExon_active_df.set_index('1st exon ID')[["CAGE", "proActiv"]] + 1)

header = ["kallisto TSS", "TPM", "strand", "proActiv ID",
         "internalPromoter", "promoterPosition", "GM12878SJ.out", "K562SJ.out"]

kallistovsproActiv_GM12878_df = pd.read_csv("/home/naoto/TSSchoiceQTL/data/kallisto_proActiv_ENCODE_GM12878.tsv",
                                    sep = "\t", names = header)

kallistovsproActiv_GM12878_df['proActiv ID'] = kallistovsproActiv_GM12878_df['proActiv ID'].astype(int).astype(str)

# kallisto
kallisto_GM12878_df = \
kallistovsproActiv_GM12878_df.groupby('proActiv ID')[["TPM"]].sum().reset_index()

# Merge
kallistovsproActiv_GM12878_plot_df = \
pd.merge(proActiv_GM12878_df, kallisto_GM12878_df, on = 'proActiv ID')

# Remove inactive promoter
kallistovsproActiv_GM12878_plot_df = \
kallistovsproActiv_GM12878_plot_df[
   ~((kallistovsproActiv_GM12878_plot_df['GM12878SJ.out'] == 0) &\
   (kallistovsproActiv_GM12878_plot_df['TPM'] == 0))
]

# Log transformation
kallistovsproActiv_GM12878_plot_df['Log GM12878SJ.out'] = \
np.log10(kallistovsproActiv_GM12878_plot_df['GM12878SJ.out'] + 1)
kallistovsproActiv_GM12878_plot_df['Log TPM'] = \
np.log10(kallistovsproActiv_GM12878_plot_df['TPM'] + 1)

# external promoters
kallistovsproActiv_GM12878_plot_external_df = \
kallistovsproActiv_GM12878_plot_df[
    kallistovsproActiv_GM12878_plot_df['internalPromoter'] == False
]

# Import
header_kallisto = ['chr', 'start', 'end', "kallisto TSS", 'kallisto score', 'strand', '1st exon ID']
kallisto_GM12878_df = pd.read_csv("/home/naoto/TSSchoiceQTL/data/kallisto_merged_sorted_1stExon.bed", sep = "\t",
                             names = header_kallisto)

 # kallisto score of 1st exon
kallisto_GM12878_1stExon_df = \
kallisto_GM12878_df.groupby(['1st exon ID', 'kallisto TSS'])[["kallisto score"]].sum().reset_index()[["1st exon ID", "kallisto score"]]

# Merge
CAGEvskallisto_GM12878_1stExon_df = pd.merge(CAGE_GM12878_1stExon_df, kallisto_GM12878_1stExon_df, on = '1st exon ID')
CAGEvskallisto_GM12878_1stExon_df.columns = ['1st exon ID', 'CAGE', 'kallisto']

# Remove inactive promoter
CAGEvskallisto_GM12878_1stExon_active_df = \
CAGEvskallisto_GM12878_1stExon_df[
    ~((CAGEvskallisto_GM12878_1stExon_df['CAGE'] == 0) & \
    (CAGEvskallisto_GM12878_1stExon_df['kallisto'] == 0))
]

CAGEvskallisto_GM12878_1stExon_active_plot_df = \
np.log10(CAGEvskallisto_GM12878_1stExon_active_df.set_index('1st exon ID')[["CAGE", "kallisto"]] + 1)

# Remove internal promoter
CAGEvskallisto_GM12878_1stExon_active_external_df = \
pd.merge(CAGEvskallisto_GM12878_1stExon_active_df, proActiv_GM12878_df[['1st exon ID', 'internalPromoter']].drop_duplicates(), on = '1st exon ID')
CAGEvskallisto_GM12878_1stExon_active_external_df = \
CAGEvskallisto_GM12878_1stExon_active_external_df[
    CAGEvskallisto_GM12878_1stExon_active_external_df['internalPromoter'] == False
]

# Log transformation
CAGEvskallisto_GM12878_1stExon_active_external_plot_df = \
np.log10(CAGEvskallisto_GM12878_1stExon_active_external_df.set_index("1st exon ID")[["CAGE", "kallisto"]] + 1)


# Merge
CAGE_kallisto_proActiv_df = \
pd.merge(CAGEvskallisto_GM12878_1stExon_active_plot_df, CAGEvsproActiv_GM12878_1stExon_active_plot_df[["proActiv"]],
        left_index = True, right_index = True)

# Drop inactive promoters
CAGE_kallisto_proActiv_df = \
CAGE_kallisto_proActiv_df[
    ~((CAGE_kallisto_proActiv_df['CAGE'] == 0) &
    (CAGE_kallisto_proActiv_df['proActiv'] == 0) &
    (CAGE_kallisto_proActiv_df['kallisto'] == 0))
]

g = sns.pairplot(
    data = CAGE_kallisto_proActiv_df,
    corner = True,
    kind = 'reg', diag_kind = 'hist',
    plot_kws = {'scatter_kws':{"alpha": 0.05, "s": 0.5}, 'color':'#455681'},
    diag_kws = {'color':'#455681'})
plt.savefig("/home/naoto/TSSchoiceQTL/png_svg/CAGE_kallisto_proActiv_AllPromoters_GM12878.png", dpi = 400)