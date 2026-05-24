#!/bin/bash
#SBATCH --partition=
#SBATCH --job-name=multivcf
#SBATCH --cpus-per-task=48
#SBATCH --mem=244G
#SBATCH --time=23:59:59

#Comando de DRAGEN para fusionar VCF
/opt/edico/bin/dragen  --ref-dir /data/hg38_graph \
     --ht-reference /data/GRCh38.fa \
     --enable-gvcf-genotyper-iterative true \
     --gvcfs-to-msvcf true \
     --variant /data/vcf_original/VCF_original_1.vcf.gz \
     --variant /data/vcf_original/VCF_original_2.vcf.gz \
     --output-directory data/output \
     --output-file-prefix msVCF

