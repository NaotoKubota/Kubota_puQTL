#!/bin/bash

# conda activate stringtie

FASTQPATH=/home/naoto/data4/GEUVADIS/fastq
STARPATH=/home/naoto/data4/GEUVADIS/STAR/mapping

# make GTF
for ID in `ls ${FASTQPATH}/fastp_json/ | grep -e HG -e NA | sed -e 's/_fastp\.json//g'`
do

	stringtie \
	-p 12 \
	--conservative \
	-G /home/naoto/STAR/gtf/Homo_sapiens.GRCh38.104.gtf \
	-o ${STARPATH}/${ID}/2nd/${ID}.stringtie.gtf \
	${STARPATH}/${ID}/2nd/${ID}.sort.bam

done

# merge GTF
mkdir -p ${STARPATH}/gtf && \
GTFLIST=`ls ${STARPATH}/*/2nd/*.stringtie.gtf` && \
stringtie --merge \
-p 24 \
-G /home/naoto/STAR/gtf/Homo_sapiens.GRCh38.104.gtf \
-o ${STARPATH}/gtf/GEUVADIS_merge.gtf \
${GTFLIST}
