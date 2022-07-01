#!/bin/bash

VCFPATH=/mnt/houman/narumi/geuvadis_RNAdata_usedata

# Liftover VCF hg19 to hg38
for ID in `ls /home/naoto/data4/GEUVADIS/fastq/fastp_json/ | grep -e HG -e NA | sed -e 's/_fastp\.json//g'`
do

	if [ ! -f /home/naoto/data4/GEUVADIS/vcf/${ID}/${ID}.hg38.vcf.gz ]; then

		echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${ID}" > header.txt && \
		mkdir -p /home/naoto/data4/GEUVADIS/vcf/${ID} && \
		cat <(grep -e "##" ${VCFPATH}/${ID}/${ID}.vcf) header.txt <(grep -v "##" ${VCFPATH}/${ID}/${ID}.vcf) \
		> /home/naoto/data4/GEUVADIS/vcf/${ID}/${ID}.hg19.vcf && \
		java -jar ~/src/picard.jar LiftoverVcf \
		-I /home/naoto/data4/GEUVADIS/vcf/${ID}/${ID}.hg19.vcf \
		-O /home/naoto/data4/GEUVADIS/vcf/${ID}/${ID}.hg38.vcf \
		-C /home/naoto/data4/GEUVADIS/vcf/hg38/b37ToHg38.over.chain \
		--REJECT /home/naoto/data4/GEUVADIS/vcf/${ID}/${ID}.hg38.unlifted.vcf \
		-R /home/naoto/data6/twins1/RefGenome/hg38.fa \
		--WARN_ON_MISSING_CONTIG true \
		--TMP_DIR ~/tmp/ && \
		bgzip /home/naoto/data4/GEUVADIS/vcf/${ID}/${ID}.hg38.vcf && \
		bgzip /home/naoto/data4/GEUVADIS/vcf/${ID}/${ID}.hg38.unlifted.vcf && \
		tabix -p vcf /home/naoto/data4/GEUVADIS/vcf/${ID}/${ID}.hg38.vcf.gz
		
	fi

done

rm -rf header.txt

# Merge vcf
VCFLIST=`ls /home/naoto/data4/GEUVADIS/vcf/???????/*.hg38.vcf.gz` && \
mkdir -p /home/naoto/data4/GEUVADIS/vcf/merge && \
bcftools merge ${VCFLIST} -o /home/naoto/data4/GEUVADIS/vcf/merge/merge.vcf.gz -O z

# MAF > 0.01
bcftools +missing2ref \
/mnt/data5/naoto/GEUVADIS/vcf/merge.vcf.gz | \
bcftools +fill-tags \
-- -t AF | \
bcftools view --exclude 'AF[0]<0.01 | AF[0]>0.99' | \
bgzip -c \
> /mnt/data5/naoto/GEUVADIS/vcf/merge.MAF001.vcf.gz && \
tabix -p vcf /mnt/data5/naoto/GEUVADIS/vcf/merge.MAF001.vcf.gz
