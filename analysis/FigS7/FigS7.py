import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')
import numpy as np

puQTL_not_eQTL_df = pd.read_csv("../data/puQTL_not_eQTL.bed", sep = "\t", names = ["chr", "start", "end", "rsid", "prmtr", "gene", "strand"])
puQTL_overlap_eQTL_df = pd.read_csv("../data/puQTL_overlap_eQTL.bed", sep = "\t", names = ["chr", "start", "end", "rsid", "prmtr", "gene", "strand"])

activity_df = pd.read_csv("../data/promoter_activity.bed.gz", sep = "\t")
activity_df = activity_df.set_index("prmtr.promoterId")

activity_score_df = activity_df.iloc[:, 5:]

activity_score_not_eQTL_df = \
activity_score_df[
    activity_score_df.index.isin(puQTL_not_eQTL_df['prmtr'])
]

activity_score_overlap_eQTL_df = \
activity_score_df[
    activity_score_df.index.isin(puQTL_overlap_eQTL_df['prmtr'])
]

activity_score_not_eQTL_df['mean'] = np.mean(activity_score_not_eQTL_df, axis = 1)
activity_score_not_eQTL_df['median'] = np.median(activity_score_not_eQTL_df, axis = 1)

activity_score_overlap_eQTL_df['mean'] = np.mean(activity_score_overlap_eQTL_df, axis = 1)
activity_score_overlap_eQTL_df['median'] = np.median(activity_score_overlap_eQTL_df, axis = 1)

activity_score_overlap_eQTL_df['group'] = "puQTL and eQTL genes"
activity_score_not_eQTL_df['group'] = "puQTL but not eQTL genes"

activity_score_df = pd.concat(
        [activity_score_overlap_eQTL_df[["group", "mean", "median"]],
        activity_score_not_eQTL_df[["group", "mean", "median"]]],
        axis = 0
)

# figure
fig = plt.figure(figsize = (4, 4))
ax = fig.add_subplot(111)

# violin plot
## median
sns.violinplot(data = activity_score_df, x = "group", y = "median", ax = ax, palette = "Set2", cut = 0)

sns.despine()
plt.tight_layout()

plt.savefig("../figure/puQTL_activitylevel_median_violin.png", dpi = 400)