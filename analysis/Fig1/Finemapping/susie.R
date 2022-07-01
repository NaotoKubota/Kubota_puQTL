## usage ##

# Rscript susie.R phenotype_ID

###########

library(susieR)
set.seed(1)
setwd("/wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/SNP-SV/susie/")
args <- commandArgs(trailingOnly = T)

# Read genotype table
genotype_table <- paste(args[1], ".genotype.txt", sep = "")
d <- read.table(genotype_table, row.names = 1)
# drop rows with the same values in all columns
keep <- apply(d, 1, function(x) length(unique(x[!is.na(x)])) != 1)
d <- d[keep, ]

# Read phenotype table
phenotype_table <- paste(args[1], ".txt", sep = "")
col.names <- list("variant", "beta", "sebeta")
p <- read.table(phenotype_table, row.names = 1, col.names = col.names)
# drop rows with zero sebeta
p <- p[p$sebeta != 0, ]

# Merge genotype and phenotype tables
dp <- merge(d, p, by=0)

# Correlation matrix
d <- d[dp$Row.names, ]
d.t <- t(d)
R <- cor(d.t)

# Calculate z scores (beta / sebeta)
z_scores <- dp$beta / dp$sebeta

# Fine-mapping by SuSiE
fitted_rss <- susie_rss(z_scores, R, L=10)

# Write a file
output_file <- paste(args[1], ".result.txt", sep = "")
write.table(fitted_rss$pip, output_file, sep = "\t", quote = F, col.names = F, append = F)
