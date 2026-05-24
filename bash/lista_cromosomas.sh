#!/bin/bash
#SBATCH --job-name=lista_chrm     # Job name
#SBATCH --output=lista_chrm_%j.out       # Output file (%j expands to job ID)
#SBATCH --error=lista_chrm_%j.err        # Error file
#SBATCH --time=02:00:00             # Walltime (3 minute)
#SBATCH --ntasks=1                  # Number of tasks (processes)
#SBATCH --cpus-per-task=4           # Number of CPU cores per task
#SBATCH --mem=8G                 # Memory per node
#SBATCH --partition=       # Queue/partition (adjust to your system)

#Script para listar todas las regiones cromosómicas de un archivo dentro de las "válidas" (autosomas, cromosomas sexuales, cromosoma M)

conda activate tfg-env
bcftools query -f '%CHROM\n' "data/msVCF_merged.vcf.gz"| sort -u | grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$'> "data/chromosomes.txt"
echo "Chromosomes list generated:"
cat "data/chromosomes.txt"
