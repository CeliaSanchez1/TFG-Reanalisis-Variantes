#!/bin/bash
#SBATCH --job-name=vcf_annot
#SBATCH --output=chr_%A_%a.out
#SBATCH --error=chr%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --array=0-23 
#SBATCH --partition=        # Queue/partition (adjust to your system)

ml Java
ml Miniconda3/24.7.1-0
conda activate python3.14

snpeff_folder="/data/snpEff"
outdir="/data/cromosomas_anotados"
input="/data/cromosomas_separados"

SNPEFF_JAR="$snpeff_folder/snpEff.jar"
SnpSift_JAR="$snpeff_folder/SnpSift.jar"
dbsnpFile=$(cat "$snpeff_folder/dbs/dbsnp.version")
dbsnp_path="$snpeff_folder/dbs/$dbsnpFile"
dbnsfpFile=$(cat "$snpeff_folder/dbs/dbnsfp.version")
dbnsfp_path=$(ls $snpeff_folder/dbs/dbNSFP*.txt.gz | head -n1)
ClinVarFile=$(cat "$snpeff_folder/dbs/clinvar.version")
ClinVar_path="$snpeff_folder/dbs/$ClinVarFile"

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "/data/chromosomes.txt")

if [ -z "$CHR" ]; then
  echo "No chromosome assigned to task ${SLURM_ARRAY_TASK_ID}"
  exit 0
fi
echo "Procesando $CHR"

conda run -n python3.14 scritps_python/anotacion_snpEff.py "$input/${CHR}.vcf.gz" "$output/snpEff/vcf_${CHR}.vcf.gz" 
java -jar "$snpeff_folder/SnpSift.jar" annotate "$dbsnp_path" "$output/snpeff/vcf_${CHR}.vcf.gz" | bgzip > "$output/dbsnp/vcf_${CHR}.vcf.gz"
tabix -p vcf "$output/cromosomas_anotados/dbsnp/vcf_${CHR}.vcf.gz"
java -jar "$snpeff_folder/SnpSift.jar" dbnsfp -v -db "$dbnsfp_path" "$output/dbsnp/vcf_${CHR}.vcf.gz" | bgzip > "$output/dbnsfp/vcf_${CHR}.vcf.gz"
tabix -p vcf "$output/dbnsfp/vcf_${CHR}.vcf.gz"
java -jar "$snpeff_folder/SnpSift.jar" annotate "$ClinVar_path" "$output/dbnsfp/vcf_${CHR}.vcf.gz" | bgzip > "$output/final/vcf_${CHR}.vcf.gz"
tabix -p vcf "$output/final/vcf_${CHR}.vcf.gz"

#Elimino las carpetas intermedias y la lista de cromosomas
rm "/data/chromosomes.txt"
rm -r "$output/snpeff"
rm -r "$output/dbsnp"
rm -r "$output/dbnsfp"
