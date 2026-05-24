#!/bin/bash
#SBATCH --job-name=fusionar_vcf     # Job name
#SBATCH --output=vcf_%j.out       # Output file (%j expands to job ID)
#SBATCH --error=vcf_%j.err        # Error file
#SBATCH --time=02:00:00             # Walltime (3 minute)
#SBATCH --ntasks=1                  # Number of tasks (processes)
#SBATCH --cpus-per-task=4           # Number of CPU cores per task
#SBATCH --mem=8G                 # Memory per node
#SBATCH --partition=        # Queue/partition (adjust to your system)

ml Java
conda activate tfg-env

SnpSift_JAR=/data/snpEff/SnpSift.jar
Mitomap_disease_db=/data/snpEff/variant_dbs/disease.vcf.gz

#Comprimir e indexar la db:
bgzip "$Mitomap_disease_db"
tabix -p vcf "$Mitomap_disease_db.gz"

#Edito el VCF porque ahí mis CHROM aparecen como chrM y en MITOMAP como MT
zcat "/data/chrM.vcf.gz" | sed 's/^chrM\t/MT\t/' | bgzip > /data/vcf_MT.vcf.gz
tabix -p vcf "/data/vcf_MT.vcf.gz"

java -jar "$SnpSift_JAR" annotate -v -db "$Mitomap_disease_db.gz" "/data/vcf_MT.vcf.gz" \
| tee "/data/chrM_final.vcf" | bgzip > "/data/chrM_final.vcf.gz"






