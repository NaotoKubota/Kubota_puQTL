import pandas as pd

TMM_df = pd.read_csv("/home/naoto/data4/GEUVADIS/STAR/featureCounts/filtered_counts_stringtie_TMM.txt",
                    sep = "\t")
TMM_df = TMM_df.rename(columns={"Unnamed: 0": "gid"})

# Make Gene bed
ENSG_GTF_df = pd.read_csv("/mnt/data4/naoto/GEUVADIS/STAR/mapping/gtf/GEUVADIS_merge.strands.gtf",
                         sep = "\t", header = None, comment = "#")

ENSG_bed_df = ENSG_GTF_df[ENSG_GTF_df[2] == "transcript"][[0, 3, 4, 6, 8]]
ENSG_bed_df['ENSG'] = ENSG_bed_df[8].str.split(";", expand = True)[0].str.split('"', expand = True)[1]
ENSG_bed_df = ENSG_bed_df[[0, 3, 4, 6, "ENSG"]]
ENSG_bed_df.columns = ["chr", "start", "end", "strand", "ENSG"]
ENSG_bed_df['chr'] = "chr" + ENSG_bed_df['chr'].astype(str)
ENSG_bed_df["ENSG_2"] = ENSG_bed_df["ENSG"]
ENSG_bed_df.columns = ["#Chr", "start", "end", "strand", "pid", "gid"]
ENSG_bed_df = ENSG_bed_df[["#Chr", "start", "end", "pid", "gid", "strand",]]
ENSG_bed_df = ENSG_bed_df.drop_duplicates(subset = "pid")

TMM_merged_df = pd.merge(ENSG_bed_df, TMM_df, on = "gid")
TMM_merged_df = TMM_merged_df.sort_values(["#Chr", "start"], ascending = [True, True])

TMM_merged_df.to_csv("/home/naoto/TSSchoiceQTL/data/eQTL/bed/TMM.bed",
                    sep = "\t", index = False)
