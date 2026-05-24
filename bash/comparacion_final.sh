#!/bin/bash
#SBATCH --job-name=comparar_vcf     # Job name
#SBATCH --output=vcf_%A_%a.out
#SBATCH --error=vcf_%A_%a.err
#SBATCH --time=02:00:00             # Walltime (3 minute)
#SBATCH --ntasks=1                  # Number of tasks (processes)
#SBATCH --cpus-per-task=4           # Number of CPU cores per task
#SBATCH --mem=8G                 # Memory per node
#SBATCH --partition=intel_skylake           # Queue/partition (adjust to your system)
#SBATCH --array=0-23   # 24 tareas

ml Miniconda3/24.7.1-0
conda activate tfg-env

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$3/chromosomes.txt")
echo "Procesando $CHR"

# Construir nombres de archivo
OLD_VCF="$1/vcf_${CHR}.vcf.gz"
NEW_VCF="$2/cromosomas_anotados/final/vcf_${CHR}.vcf.gz"

wget -O data/phenotype_to_genes.txt http://purl.obolibrary.org/obo/hp/phenotype_to_genes.txt
conda run -n tfg-env python scritps_python/comparar.py "data/phenotype_to_genes.txt" "$OLD_VCF" "$NEW_VCF" "$2" "$2/Union_HPO/${CHR}.tsv" "$2/Union_HPO_reducido/${CHR}.tsv"



