import pandas as pd

puQTL_df = pd.read_csv("/home/naoto/TSSchoiceQTL/data/puQTL/bed/proActiv_GEUVADIS_merged_stringtie_QTLtools.bed.gz", sep = "\t")
puQTL_score_df = puQTL_df.T.iloc[6:].astype(float)
puQTL_score_describe_df = puQTL_score_df.describe()

# Remove the promoters of which the 4th quantile of score is 0
puQTL_score_expressedTSS_l = \
list(puQTL_score_describe_df.T[puQTL_score_describe_df.T['75%'] != 0].index)
puQTL_expressed_df = puQTL_df.T[puQTL_score_expressedTSS_l].T

# Export
puQTL_expressed_df.to_csv("/home/naoto/TSSchoiceQTL/data/puQTL/bed/promoter_activity.bed", sep = "\t", header = True, index = False)
