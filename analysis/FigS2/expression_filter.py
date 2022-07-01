import pandas as pd

# +
count_df = pd.read_csv("/home/naoto/data4/GEUVADIS/STAR/featureCounts/all_counts_stringtie.txt",
                    sep = "\t", comment = "#", index_col = 0)

ID = list(pd.DataFrame(count_df.columns[5:])[0].str.split("/", expand = True)[7])
col = list(count_df.columns[0:5]) + ID
count_df.columns = col
# -

count_simple_df = count_df.iloc[:, 5:]
count_simple_df = count_simple_df.T[count_simple_df.sum() != 0].T

# TPM

# +
gene_length = count_df['Length']
def normalize_per_kilobase(df, gene_length):
    df_tmp = df.copy()
    df_tmp = (df_tmp.T / gene_length * 10 **3).T
    return df_tmp

def normalize_per_million_reads(df):
    sum_count = df.sum()
    return 10**6 * df / sum_count

# TPM
def normalize_tpm(df, gene_length):
    df_tmp = df.copy()
    df_tmp = normalize_per_kilobase(df_tmp, gene_length)
    df_tmp = normalize_per_million_reads(df_tmp)
    return df_tmp


# -

count_tpm_df = normalize_tpm(count_simple_df, gene_length)

# >0.1 TPM in at least 20% of samples
TPM_bool = (count_tpm_df > 0.1).T.sum() > 438*0.8
# ≥6 reads in at least 20% of samples
Reads_bool = (count_simple_df > 6).T.sum() > 438*0.8

count_filtered_df = count_simple_df[TPM_bool]
count_filtered_df = count_filtered_df[Reads_bool]

# Drop mitochondria genes
count_filtered_df = count_filtered_df[~count_filtered_df.index.str.contains("MT")]

count_filtered_df.to_csv("/home/naoto/data4/GEUVADIS/STAR/featureCounts/filtered_counts_stringtie.txt",
                        sep = "\t", index = True, header = True)
