library(proActiv)

## From GTF file path
gtf.file <- '/wakasugi/mnt/data4/naoto/GEUVADIS/STAR/mapping/gtf/GEUVADIS_merge.strands.gtf'
promoterAnnotation.geuvadis.merge.subset <- preparePromoterAnnotation(file = gtf.file, species = 'Homo_sapiens')
## List of STAR junction files as input
files <- list.files('/wakasugi/mnt/data4/naoto/GEUVADIS/STAR/junctions', full.names=T)
## Vector describing experimental condition
condition <- rep(c('LCL'), each=1)
## Promoter annotation for human genome GENCODE v104
promoterAnnotation <- promoterAnnotation.geuvadis.merge.subset

result <- proActiv(files = files, promoterAnnotation = promoterAnnotation)

## Removes single-exon transcripts / promoters by eliminating promoter counts that are NA 
result <- result[complete.cases(assays(result)$promoterCounts),]

## export result
result_df <- rowData(result)
write.table(result_df, "/wakasugi/home/naoto/TSSchoiceQTL/data/proActiv/result_rowData_GEUVADIS.tsv", sep = "\t", quote = F)
result_absolutePromoterActivity_df <- assays(result)$absolutePromoterActivity
write.table(result_absolutePromoterActivity_df, "/wakasugi/home/naoto/TSSchoiceQTL/data/proActiv/result_absolutePromoterActivity_GEUVADIS.tsv", sep = "\t", quote = F)
