#!/bin/bash
#SBATCH --job-name=fusionar_vcf     # Job name
#SBATCH --output=vcf_%j.out       # Output file (%j expands to job ID)
#SBATCH --error=vcf_%j.err        # Error file
#SBATCH --time=02:00:00             # Walltime (3 minute)
#SBATCH --ntasks=1                  # Number of tasks (processes)
#SBATCH --cpus-per-task=4           # Number of CPU cores per task
#SBATCH --mem=8G                 # Memory per node
#SBATCH --partition=         # Queue/partition (adjust to your system)

#Wrapper para comparar versiones de la anotación empleando la función de python
conda activate tfg-env

OLD_VCF="/data/Comparacion_versiones/anotacion-version-anterior/vcf_anotado.vcf.gz"
NEW_VCF="/data/Comparacion_versiones/anotacion-ultima-version/vcf_anotado.vcf.gz"
OUTDIR="/data/Comparacion_versiones/resultados.gz"
mkdir -p "$OUTDIR"

# COMPARACIÓN
conda run -n tfg-env python comparar_versiones_anotacion.py "$OLD_VCF" "$NEW_VCF" "$OUTDIR" 

