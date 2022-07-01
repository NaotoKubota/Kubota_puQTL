#!/bin/bash

# conda activate star

FASTQPATH=/home/naoto/data4/GEUVADIS/fastq
STARPATH=/home/naoto/data4/GEUVADIS/STAR/mapping

### STAR 2-pass mapping ###
# 1st mapping
for ID in `ls ${FASTQPATH}/fastp_json/ | grep -e HG -e NA | sed -e 's/_fastp\.json//g'`
do

    if [ ! -f ${STARPATH}/${ID}/1st/${ID}_SJ.out.tab ]; then
    
        mkdir -p ${STARPATH}/${ID}/1st/ && \
        STAR \
        --runThreadN 32 \
        --genomeDir /home/naoto/data4/GEUVADIS/STAR/index_hg38_ensembl_v104 \
        --readFilesIn \
        ${FASTQPATH}/${ID}/${ID}_1_filtered.fastq.gz \
        ${FASTQPATH}/${ID}/${ID}_2_filtered.fastq.gz \
        --readFilesCommand zcat \
        --outFilterMultimapNmax 1 \
        --outFileNamePrefix ${STARPATH}/${ID}/1st/${ID}_
    
    fi

done

# 2nd mapping
SJLIST=`ls ${STARPATH}/*/1st/*tab`
for ID in `ls ${FASTQPATH} | grep -e HG -e NA | sed -e 's/\///g'`
do

    if [ ! -f ${STARPATH}/${ID}/2nd/${ID}_SJ.out.tab ]; then

        mkdir -p ${STARPATH}/${ID}/2nd/ && \
        STAR \
        --runThreadN 32 \
        --genomeDir /home/naoto/data4/GEUVADIS/STAR/index_hg38_ensembl_v104 \
        --readFilesIn \
        ${FASTQPATH}/${ID}/${ID}_1_filtered.fastq.gz \
        ${FASTQPATH}/${ID}/${ID}_2_filtered.fastq.gz \
        --readFilesCommand zcat \
        --outFilterMultimapNmax 1 \
        --sjdbFileChrStartEnd ${SJLIST} \
        --limitSjdbInsertNsj 10000000 \
        --outFileNamePrefix ${STARPATH}/${ID}/2nd/${ID}_ && \
        samtools sort \
        -@ 16 -O bam \
        -o ${STARPATH}/${ID}/2nd/${ID}.sort.bam \
        ${STARPATH}/${ID}/2nd/${ID}_Aligned.out.sam && \
        samtools index ${STARPATH}/${ID}/2nd/${ID}.sort.bam

    fi

done
