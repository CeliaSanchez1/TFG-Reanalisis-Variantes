#!/bin/bash
#SBATCH --job-name=workflow    # Job name
#SBATCH --output=workflow_%j.out       # Output file (%j expands to job ID)
#SBATCH --error=workflow_%j.err        # Error file
#SBATCH --time=10:00:00             # Walltime (3 minute)
#SBATCH --ntasks=1                  # Number of tasks (processes)
#SBATCH --cpus-per-task=4           # Number of CPU cores per task
#SBATCH --mem=8G                 # Memory per node
#SBATCH --partition=         # Queue/partition (adjust to your system)

conda activate tfg-env
set -euo pipefail

#Directorios
OUTDIR= #Ruta al directorio donde se guardarán los archivos  generados en este workflow
ref_dir= #Ruta al directorio del genoma de referencia
ht_ref= #Ruta al archivo de referencia para el genotipado iterativo
old_comparacion_dir= #Ruta a la carpeta de los tsv viejos que contienen las muestras modificadas de una anotación a otra-->ESTO ES CARPETA, CON 23 CHROMS
old_comparacion_CHRM_dir= #Ruta a la carpeta de los tsv viejos que contienen las muestras del CHRM modificadas de una anotación a otra-->ESTO ES RUTA A UN ÚNICO VCF
batch_dir= #Ruta a la carpeta del nuevo lote de VCFs
snpEff_dir= #Ruta a la carpeta que contiene ejecutables y dbs de snpEff
db_dir= #Ruta a la carpeta con las bases de datos (normalmente dentro de la de snpEff)
msVCF_dir= #Ruta al msVCF actual

mkdir -p "$OUTDIR"

##CADA VEZ QUE SE RECIBE UN LOTE DE MUESTRAS SE EJECUTA LO SIGUIENTE:

#1-Unir las muestras al msVCF que ya tengo
sbatch iterar_msVCF_lista.sh "$batch_dir" "$ref_dir" "$ht_ref" "$msVCF_dir" "$OUTDIR" 

#2-Comprobar que las versiones de las db sean las más actualizadas:
sbatch descarga_db.sh "$snpEff_dir" "$db_dir" #AQUÍ BUSCAR SACAR UNA FLAG PARA LA REANOTACIÓN

#3-Buscar las regiones cromosómicas en mi msVCF
sbatch lista_cromosomas.sh "$OUTDIR"

#4-Separar en 25 msVCF, uno por cromosoma
sbatch separar_cromosomas.sh "$OUTDIR"

#5-Hacer la reanotación, borro el cromosoma M de la lista y lo gestiono por separado
sed -i '/chrM/d' "$OUTDIR/chromosomes.txt"
sbatch anotar_cromosomas_paralelo.sh "$snpEff_dir" "$db_dir" "$OUTDIR" 
sabtch anotar_chrM.sh "$snpEff_dir" "$db_dir" "$OUTDIR" 

#6-Comparar con versiones anteriores de la anotación y añadir info de los HPO de los pacientes
sbatch comparacion.sh "$old_comparacion_dir" "$OUTDIR" 
sbatch comparacion_chrM.sh "$old_comparacion_CHRM_dir" "$OUTDIR" 

#7-Volver a unir los archivos de los distintos cromosomas entre sí para obtener un tsv final
sbatch fusion_columnas.sh "$OUTDIR"
