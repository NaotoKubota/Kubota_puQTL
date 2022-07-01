####################### Must run in docker container! #########################
# docker run -v /:/wakasugi --rm -it naotokubota/coloc-locuscomparer:1.0 bash
###############################################################################

# import library
library("coloc")

# set args
args <- commandArgs(trailingOnly = T)
input_1 <- args[1]
input_2 <- args[2]
PROPOTION <- as.numeric(args[3])
output <- args[4]

# read two input tables
d1 <- read.table(file=input_1, header=T, as.is=T) # puQTL
d2 <- read.table(file=input_2, header=T, as.is=T) # GWAS (Case-Control)
# read MAF table
MAF <- read.table(file="/wakasugi/mnt/data5/naoto/GEUVADIS/vcf/MAF_table.txt", header=T, as.is=T)

# merge
d1_MAF <- merge(d1, MAF, by="rs_id", all=FALSE)
input <- merge(d1_MAF, d2, by="rs_id", all=FALSE, suffixes=c("_puQTL", "_gwas"))
input <- na.omit(input)
input <- transform(input, MAF = as.numeric(MAF))
input <- na.omit(input)
input <- transform(input, MAF = as.numeric(MAF))

# run coloc
result <- coloc.abf(
dataset1=list(pvalues=input$pval_nominal, type="quant", N=nrow(d1)),
dataset2=list(pvalues=input$pval, type="cc", s=PROPOTION, N=nrow(d2)),
MAF=input$MAF)
df <- result[[1]]
df <- as.data.frame(t(data.frame(df)))
df <- transform(df, promoter=c(input$prmtr_id[1]))
df <- transform(df, gwas=c(input$gwas[1]))

# save
write.table(df, output, sep = "\t", quote = F, col.names=F, row.names=F, append = T)