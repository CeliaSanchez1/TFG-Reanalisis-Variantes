#!/bin/bash
#SBATCH --partition=
#SBATCH --job-name=multivcf
#SBATCH --cpus-per-task=48
#SBATCH --mem=244G
#SBATCH --time=23:59:59

#Todas las muestras almacenadas en la carpeta VCF_originales.
#Listo las muestras en un txt para evitar pasar tantos argumentos en el comando de DRAGEN.
find "/data/VCF_originales" -type f -name "*.vcf.gz" | sort > "data/lista_muestras_vcf.txt"

/opt/edico/bin/dragen  --ref-dir /data/hg38_graph \
     --ht-reference /data/GRCh38.fa \
     --enable-gvcf-genotyper-iterative true \
     --gvcfs-to-msvcf true \
     --variant-list "data/lista_muestras_vcf.txt" \
     --output-directory "data/fusion.vcf.gz" \
     --output-file-prefix msVCF

rm "data/lista_muestras_vcf.txt" #Borro la lista intermedia de muestras una vez he realizado la fusión
