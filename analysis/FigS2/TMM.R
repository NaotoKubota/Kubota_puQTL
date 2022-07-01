### in docker ###
# docker run --rm -v /:/wakasugi broadinstitute/gtex_eqtl:V8 bash

# Import
library("edgeR")

count <- read.table("/wakasugi/mnt/data4/naoto/GEUVADIS/STAR/featureCounts/filtered_counts.txt", sep = "\t", header = T, row.names = 1)
count <- as.matrix(count)

# TMM Normalization
d <- DGEList(counts = count)
d <- calcNormFactors(d, method = "TMM")
Norm.Factors <- d$samples$norm.factors
Lib.Size <- d$samples$lib.size
count.norm <- t(apply(count, 1, function(x) x / (Norm.Factors * Lib.Size) * 10^6))

# Save
write.table(count.norm, "/wakasugi/mnt/data4/naoto/GEUVADIS/STAR/featureCounts/filtered_counts_TMM.txt", sep = "\t", quote = F, col.names = NA)
