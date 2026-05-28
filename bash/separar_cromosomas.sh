#!/bin/bash
#SBATCH --job-name=separar_cromosomas     # Job name
#SBATCH --output=separar_cromosomas_%j.out       # Output file (%j expands to job ID)
#SBATCH --error=separar_cromosomas_%j.err        # Error file
#SBATCH --time=2:00:00             # Walltime 
#SBATCH --array=1-25          # Number of tasks (processes)
#SBATCH --cpus-per-task=1           # Number of CPU cores per task
#SBATCH --mem=4G                 # Memory per node
#SBATCH --partition=           # Queue/partition 

#Script para dividir el msVCF en otros msVCF a partir de la lista de regiones cromosómicas generada previamente

conda activate tfg-env

mkdir -p "data/cromosomas_separados"
CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$data/chromosomes.txt")
if [ -z "$CHR" ]; then
  echo "No chromosome assigned to task ${SLURM_ARRAY_TASK_ID}"
  exit 0
fi

echo "Processing $CHR"
bcftools view -r "$CHR" -Oz -o "data/cromosomas_separados/${CHR}.vcf.gz" "data/msVCF_merged.vcf.gz"
bcftools index -f -t "data/cromosomas_separados/${CHR}.vcf.gz"
