## docker run -v /:/wakasugi -it --name gtex_eqtl broadinstitute/gtex_eqtl:V8 bash

## usage ##

# Rscript peer.R /path/to/input /path/to/output

###########

# import library
library(peer)

# set args
args <- commandArgs(trailingOnly = T)
input_name <- args[1]
output_name <- args[2]

# import phenotype table
expr <- read.table(file = input_name, sep = "\t", header=FALSE)
expr <- as.data.frame(t(expr))

# run peer
model <- PEER()
PEER_setPhenoMean(model, as.matrix(expr))
PEER_setNk(model, 50)
PEER_getNk(model)
PEER_setNmax_iterations(model, 10000000)
PEER_update(model)

# get peer factors
factors <- PEER_getX(model)
factors_df <- data.frame(factors)
factors_df <- as.data.frame(t(factors_df))

# save
write.table(factors_df, output_name, sep = "\t", quote = F, col.names=F, append = F)
