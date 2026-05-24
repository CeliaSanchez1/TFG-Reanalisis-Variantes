#!/bin/bash
#SBATCH --job-name=comparar_vcf     # Job name
#SBATCH --output=vcf_%A_%a.out
#SBATCH --error=vcf_%A_%a.err
#SBATCH --time=23:59:59             # Walltime 
#SBATCH --ntasks=1                  # Number of tasks (processes)
#SBATCH --cpus-per-task=1          # Number of CPU cores per task
#SBATCH --mem=4G                 # Memory per node
#SBATCH --partition=        # Queue/partition (adjust to your system)

#MEJORAS RESPECTO A LA VERSIÓN ANTERIOR:
#Renombrar chrM → MT también en el contig
#Filtrar alelos simbólicos <NON_REF> para evitar NullPointerException
#Filtrar muestras sin alelo alternativo real 
#tabix con -p vcf explícito en todos los pasos

set -euo pipefail

SNPEFF_JAR="/data/.../snpEff.jar"
SNPSIFT_JAR="/.../SnpSift.jar"
MITOMAP_DB="/data/.../disease.vcf.gz"
SNPEFF_GENOME="GRCh38.99"
SNPEFF_MEM="4g"
INPUT="/data/.../chrM.vcf.gz"
WORKDIR="/data/.../anotar_chrM"
OUTPUT="${WORKDIR}/chrM_final.vcf.gz"
TMP_RENAMED="${WORKDIR}/tmp_01_renamed.vcf.gz"
TMP_CLEAN="${WORKDIR}/tmp_02_clean.vcf.gz"
TMP_SNPEFF_VCF="${WORKDIR}/tmp_03_snpeff.vcf"
TMP_SNPEFF_GZ="${WORKDIR}/tmp_03_snpeff.vcf.gz"

#Renombrar chrM → MT
CHR_MAP="${WORKDIR}/chr_map.txt"
printf "chrM\tMT\n" > "$CHR_MAP"

#Anotar
bcftools annotate --rename-chrs "$CHR_MAP" -Oz -o "$TMP_RENAMED" "$INPUT"
tabix -p vcf "$TMP_RENAMED"

#Filtrar entradas no anotables
bcftools view -e 'ALT="<NON_REF>"'"$TMP_RENAMED" | bcftools view -c 1 -Oz -o "$TMP_CLEAN"
tabix -p vcf "$TMP_CLEAN"
java -Xmx"${SNPEFF_MEM}" -jar "$SNPEFF_JAR" ann -v -nodownload "$SNPEFF_GENOME" "$TMP_CLEAN" > "$TMP_SNPEFF_VCF"

#Comprimir con bgzip (BGZF) e indexar
bgzip -f "$TMP_SNPEFF_VCF"                 
tabix -p vcf "$TMP_SNPEFF_GZ"

#Anotación clínica con SnpSift + MITOMAP
java -Xmx"${SNPEFF_MEM}" -jar "$SNPSIFT_JAR" annotate -v -tabix "$MITOMAP_DB" "$TMP_SNPEFF_GZ" | bgzip > "$OUTPUT"
tabix -p vcf "$OUTPUT"

