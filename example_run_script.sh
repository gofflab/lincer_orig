#!/usr/bin/env bash

#Mouse cortical neuron subtype only catalog build

mask_gtf=/n/rinn_data1/indexes/mouse/mm9/annotations/dk/mask.gtf
pseudogene_gtf=/n/rinn_data1/indexes/mouse/mm9/annotations/dk/pseudogenes.gtf
lincs_gtf=/n/rinn_data1/indexes/mouse/mm9/annotations/dk/lncrna.gtf

python /n/rinn_data1/users/lgoff/sw/lincer/build_transript_db.py --run-cufflinks 1 -g mm9 --ref-sequence /n/rinn_data1/indexes/igenomes/Mus_musculus/UCSC/mm9/Sequence/Chromosomes/mm9.fa \
-q unrestricted_parallel --cufflinks-queue local --sample-reads primary_asm_bams.txt --min-linc-length 200 --mask-gtf $mask_gtf --lsf-nodes 400 -p 8 --min-transcript-cov 1 \
--csf-score-threshold 100 --pseudogene $pseudogene_gtf --run-meta-assembly --external-linc-gtf $lincs_gtf primary_asms.txt > build_Arlotta.out
