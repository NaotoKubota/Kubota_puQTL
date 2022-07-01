#!/bin/bash

datapath=/mnt/houman/narumi/geuvadis_RNAdata_usedata
outputpath=/home/naoto/data4/GEUVADIS/fastq

# fastp
for sample in `ls /mnt/houman/narumi/geuvadis_RNAdata_usedata`
do
	ID=`echo ${sample} | sed -e 's/\///g'`

	if [ -f ${datapath}/${ID}/*1.fastq.gz -a -f ${datapath}/${ID}/*2.fastq.gz ]; then

		mkdir -p ${outputpath}/${ID}/

		# fastp
		fastp \
		-i ${datapath}/${ID}/*_1.fastq.gz \
		-I ${datapath}/${ID}/*_2.fastq.gz \
		-3 \
		-o ${outputpath}/${ID}/${ID}_1_filtered.fastq.gz \
		-O ${outputpath}/${ID}/${ID}_2_filtered.fastq.gz \
		-h ${outputpath}/${ID}/${ID}_fastp.html \
		-j ${outputpath}/${ID}/${ID}_fastp.json \
		-q 30 -w 16

	elif [ -f ${datapath}/${ID}/*1.fastq -a -f ${datapath}/${ID}/*2.fastq ]; then

		mkdir -p ${outputpath}/${ID}/

		# fastp
		fastp \
		-i ${datapath}/${ID}/*_1.fastq \
		-I ${datapath}/${ID}/*_2.fastq \
		-3 \
		-o ${outputpath}/${ID}/${ID}_1_filtered.fastq \
		-O ${outputpath}/${ID}/${ID}_2_filtered.fastq \
		-h ${outputpath}/${ID}/${ID}_fastp.html \
		-j ${outputpath}/${ID}/${ID}_fastp.json \
		-q 30 -w 16 && \
		gzip -c ${outputpath}/${ID}/${ID}_1_filtered.fastq > ${outputpath}/${ID}/${ID}_1_filtered.fastq.gz && \
		gzip -c ${outputpath}/${ID}/${ID}_2_filtered.fastq > ${outputpath}/${ID}/${ID}_2_filtered.fastq.gz
		
	fi

done

# multiqc
mkdir -p ${outputpath}/fastp_json/ && \
cp ${outputpath}/*/*_fastp.json ${outputpath}/fastp_json/ && \
cd ${outputpath}/fastp_json/ && \
multiqc .
