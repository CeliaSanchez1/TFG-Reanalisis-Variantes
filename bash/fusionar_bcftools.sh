#!/bin/bash
#SBATCH --job-name=fusionar_vcf     # Job name
#SBATCH --output=vcf_%j.out       # Output file (%j expands to job ID)
#SBATCH --error=vcf_%j.err        # Error file
#SBATCH --time=00:03:00             # Walltime (3 minute)
#SBATCH --ntasks=1                  # Number of tasks (processes)
#SBATCH --cpus-per-task=1           # Number of CPU cores per task
#SBATCH --mem=8G                 # Memory per node
#SBATCH --partition=        # Queue/partition

#Wrapper para fusionar archivos VCF con funciones de python
conda activate tfg-env

#Indexar archivos antes de fusionarlos 
conda run -n tfg-env bcftools index -f /data/vcf_original/VCF_original_1.vcf.gz
conda run -n tfg-env bcftools index -f /data/vcf_original/VCF_original_2.vcf.gz

#Fusionar 
conda run -n tfg-env bcftools merge /data/vcf_original/VCF_original_1.vcf.gz /data/vcf_original/VCF_original_2.vcf.gz \
 -Oz -o /data/pasos_intermedios/vcf_fusionado.vcf.gz

