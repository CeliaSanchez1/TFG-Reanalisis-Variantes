#!/bin/bash
#SBATCH --job-name=fusionar_vcf     # Job name
#SBATCH --output=vcf_%j.out       # Output file (%j expands to job ID)
#SBATCH --error=vcf_%j.err        # Error file
#SBATCH --time=10:00:00             # Walltime (ESTO PONER MÍNIMO A 5h PARA QUE DE TIEMPO A TODO EN TOTAL)
#SBATCH --ntasks=1                  # Number of tasks (processes)
#SBATCH --cpus-per-task=4           # Number of CPU cores per task
#SBATCH --mem=10G                 # Memory per node
#SBATCH --partition=       # Queue/partition (adjust to your system)

ml Java
#Script para anotar ficheros VCF aplicando las bases de datos secuencialmente

# Anotar usando dbsnp
java -jar /data/SnpSift.jar annotate \
    -v -db /data/00-All.vcf.gz \
    /data/vcf_original/vcf_original_1.vcf.gz \
    >/dat/vcf_anotado/vcf_dbsnp.vcf

#Anotar usando dbNSFP
java -jar /data/SnpSift.jar annotate dbnsfp \
    -v -db /data/dbNSFP4.1a.txt.gz \
    /data/vcf_anotado/vcf_dbsnp.vcf \
    > /data/vcf_anotado/vcf_dbsnp_dbnsfp.vcf
 
# #Anotar usando ClinVar
java -jar /data/SnpSift.jar annotate \
    -v /data/clinvar_20250209.vcf.gz \
    /data/vcf_anotado/vcf_dbsnp_dbnsfp.vcf \
    > /data/vcf_anotado/vcf_dbsnp_dbnsfp_clinvar.vcf 

#Anotar usando gnomAD
java -Xmx1g -jar /data/SnpSift.jar annotate -info AF -name gnomAD_ \
    /data/grch38.vcf.bgz \
    /data/vcf_anotado/vcf_dbsnp_dbnsfp_clinvar.vcf \
    > /data/vcf_dbsnp_dbnsfp_clinvar_gnomAD.vcf
