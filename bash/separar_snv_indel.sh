#!/bin/bash
#SBATCH --job-name=fusionar_vcf     # Job name
#SBATCH --output=vcf_%j.out       # Output file (%j expands to job ID)
#SBATCH --error=vcf_%j.err        # Error file
#SBATCH --time=00:03:00             # Walltime (3 minute)
#SBATCH --ntasks=1                  # Number of tasks (processes)
#SBATCH --cpus-per-task=1           # Number of CPU cores per task
#SBATCH --mem=8G                 # Memory per node
#SBATCH --partition=        # Queue/partition

ml Miniconda3/24.7.1-0

#Guardar solo los VCF con SNVs e indels en una versión comprimida 
conda run -n python3.14 bcftools view -v snps,indels -Oz -o /data/pasos_intermedios/VCF_SNVS.vcf.gz /data/pasos_intermedios/vcf_fusionado.vcf.gz

#Lo mismo pero en una versión sin comprimir para visualizarlo mejor
conda run -n python3.14 bcftools view -v snps,indels -Ov -o /data/pasos_intermedios/VCF_SNVS.vcf /data/pasos_intermedios/vcf_fusionado.vcf.gz
