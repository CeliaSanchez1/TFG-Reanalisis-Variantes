#!/bin/bash
#SBATCH --job-name=vcf_prueba       # Job name
#SBATCH --output=prueba_%j.out       # Output file (%j expands to job ID)
#SBATCH --error=prueba_%j.err        # Error file
#SBATCH --time=00:03:00             # Walltime (3 minute)
#SBATCH --ntasks=1                  # Number of tasks (processes)
#SBATCH --cpus-per-task=1           # Number of CPU cores per task
#SBATCH --mem=8G                 # Memory per node
#SBATCH --partition         # Queue/partition 

#Wrapper para fusionar VCF con funciones de python
conda activate tfg-env
export _JAVA_OPTIONS="-Xmx6G"   #Para poder cargar bien la base de datos al utilizar SnpEff

conda run -n tfg-env python fusionar_python.py #Esto es todo mi flujo de trabajo de pre-procesar VCF, separar SNV y SV, y anotar SNV con SnpEff. 

