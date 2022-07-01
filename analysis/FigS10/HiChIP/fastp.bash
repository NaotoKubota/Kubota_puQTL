#!/bin/bash


# fastp
for SAMPLE in SRR5831489 SRR5831490;
do
    
    datapath=/mnt/data4/naoto/puQTL/HiChIP

    # fastp
    fastp \
    -i ${datapath}/${SAMPLE}*_1.fastq.gz \
    -I ${datapath}/${SAMPLE}*_2.fastq.gz \
    -3 \
    -o ${datapath}/${SAMPLE}_R1_filtered.fastq.gz \
    -O ${datapath}/${SAMPLE}_R2_filtered.fastq.gz \
    -h ${datapath}/${SAMPLE}_fastp.html \
    -j ${datapath}/${SAMPLE}_fastp.json \
    -q 30 -w 16

done

# multiqc
docker run --rm -tv /mnt/data4/naoto/puQTL/HiChIP:/mnt/data4/naoto/puQTL/HiChIP -w /mnt/data4/naoto/puQTL/HiChIP ewels/multiqc

# Move files
mkdir -p /mnt/data4/naoto/puQTL/HiChIP/fastq/SRR5831489 && \
mkdir -p /mnt/data4/naoto/puQTL/HiChIP/fastq/SRR5831490 && \
mv \
/mnt/data4/naoto/puQTL/HiChIP/SRR5831489_R1_filtered.fastq.gz \
/mnt/data4/naoto/puQTL/HiChIP/SRR5831489_R2_filtered.fastq.gz \
/mnt/data4/naoto/puQTL/HiChIP/fastq/SRR5831489/ && \
mv \
/mnt/data4/naoto/puQTL/HiChIP/SRR5831490_R1_filtered.fastq.gz \
/mnt/data4/naoto/puQTL/HiChIP/SRR5831490_R2_filtered.fastq.gz \
/mnt/data4/naoto/puQTL/HiChIP/fastq/SRR5831490/