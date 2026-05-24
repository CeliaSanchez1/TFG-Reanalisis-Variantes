#!/bin/bash
#SBATCH --partition=
#SBATCH --job-name=multivcf
#SBATCH --cpus-per-task=48
#SBATCH --mem=244G
#SBATCH --time=23:59:59

#Script para añadir el nuevo lote al msVCF

#Listar las muestras del nuevo lote en un txt para fusionarlas después
find "/data/VCF_lote" -type f -name "*.vcf.gz" | sort > "data/lista_muestras_lote_vcf.txt"

#Fusionar el lote en un msVCF
/opt/edico/bin/dragen  --ref-dir /data/hg38_graph \
    --ht-reference /data/GRCh38.fa \
    --enable-gvcf-genotyper-iterative true --gvcfs-to-msvcf true \
    --variant-list "data/lista_muestras_lote_vcf.txt" \
    --output-directory "data/fusion_lote.vcf.gz" --output-file-prefix msVCF

#Fusionar este msVCF con el que ya tenía con todas las muestras 
/opt/edico/bin/dragen --ref-dir /data/hg38_graph \
    --ht-reference /data/GRCh38.fa \
    --enable-gvcf-genotyper-iterative true --merge-batches true \
    --variant "data/msVCF_fusion_lote.vcf.gz" \
    --variant "data/msVCF_fusion.vcf.gz" \
    --output-directory "data/merged_final.vcf.gz" \
    --output-file-prefix msVCF

rm "data/lista_muestras_lote_vcf.txt" #Borrar la lista intermedia de muestras
rm "data/msVCF_fusion_lote.vcf.gz" #Borrar el msVCF del nuevo lote una vez lo he fusionado con el otro
rm "data/msVCF_fusion.vcf.gz" #Borrar el antiguo msVCF después de añadirle todas las muestras


