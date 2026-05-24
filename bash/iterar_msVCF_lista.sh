#!/bin/bash
#SBATCH --partition=
#SBATCH --job-name=multivcf
#SBATCH --cpus-per-task=48
#SBATCH --mem=244G
#SBATCH --time=23:59:59

#Recibo como input la carpeta con el nuevo lote, las referencias del genoma y la carpeta de salida
#Guardo en la carpeta de salida el msVCF actualizado con todas las muestras
#Borro el msVCF histórico anterior

#Al recibir un nuevo lote lo primero que hago es listar sus muestras en un txt para fusionarlas después
find "/data/VCF_lote" -type f -name "*.vcf.gz" | sort > "data/lista_muestras_lote_vcf.txt"

#Luego fusiono todo el lote en un msVCF
/opt/edico/bin/dragen  --ref-dir /data/hg38_graph \
    --ht-reference /data/GRCh38.fa \
    --enable-gvcf-genotyper-iterative true --gvcfs-to-msvcf true \
    --variant-list "data/lista_muestras_lote_vcf.txt" \
    --output-directory "data/fusion_lote.vcf.gz" --output-file-prefix msVCF

#Luego fusiono este msVCF con el msVCF que ya tenía con todas las muestras 
/opt/edico/bin/dragen --ref-dir /data/hg38_graph \
    --ht-reference /data/GRCh38.fa \
    --enable-gvcf-genotyper-iterative true --merge-batches true \
    --variant "data/msVCF_fusion_lote.vcf.gz" \
    --variant "data/msVCF_fusion.vcf.gz" \
    --output-directory "data/merged_final.vcf.gz" \
    --output-file-prefix msVCF

rm "data/lista_muestras_lote_vcf.txt" #Borro la lista intermedia de muestras
rm "data/msVCF_fusion_lote.vcf.gz" #Borro el msVCF del nuevo lote una vez lo he fusionado con el otro
rm "data/msVCF_fusion.vcf.gz" #Borro el antiguo msVCF después de añadirle todas las muestras


