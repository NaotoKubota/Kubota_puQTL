#!/bin/bash

zcat /home/naoto/TSSchoiceQTL/data/puQTL/bed/promoter_activity.bed.gz | \
cut -f 7- | sed -e '1d' \
> /home/naoto/TSSchoiceQTL/data/puQTL/bed/promoter_activity_forpeer.tsv

docker run --rm -v /:/wakasugi broadinstitute/gtex_eqtl:V8 bash -c "\
Rscript /wakasugi/home/naoto/TSSchoiceQTL/scripts/peer.R \
/wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/bed/promoter_activity_forpeer.tsv \
/wakasugi/home/naoto/TSSchoiceQTL/data/puQTL/txt/promoter_activity_peerfactors.txt\
"
