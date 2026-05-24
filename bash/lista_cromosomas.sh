#!/bin/bash
#SBATCH --job-name=fusionar_vcf     # Job name
#SBATCH --output=vcf_%j.out       # Output file (%j expands to job ID)
#SBATCH --error=vcf_%j.err        # Error file
#SBATCH --time=02:00:00             # Walltime (3 minute)
#SBATCH --ntasks=1                  # Number of tasks (processes)
#SBATCH --cpus-per-task=4           # Number of CPU cores per task
#SBATCH --mem=8G                 # Memory per node
#SBATCH --partition=       # Queue/partition (adjust to your system)

#Trabajo sobre el msVCF completo, obtengo un txt que lista las regiones cromosómicas
ml Miniconda3/24.7.1-0
conda activate bcftools-env
bcftools query -f '%CHROM\n' "data/msVCF_merged.vcf.gz"| sort -u | grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$'> "data/chromosomes.txt"
echo "Chromosomes list generated:"
cat "data/chromosomes.txt"