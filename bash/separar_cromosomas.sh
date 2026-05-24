#!/bin/bash
#SBATCH --job-name=separar_cromosomas     # Job name
#SBATCH --output=vcf_%j.out       # Output file (%j expands to job ID)
#SBATCH --error=vcf_%j.err        # Error file
#SBATCH --time=2:00:00             # Walltime 
#SBATCH --array=1-25          # Number of tasks (processes)
#SBATCH --cpus-per-task=1           # Number of CPU cores per task
#SBATCH --mem=4G                 # Memory per node
#SBATCH --partition=e           # Queue/partition 

ml Miniconda3/24.7.1-0
conda activate bcftools-env
#Recibo como input el msVCF y utilizo la lista de cromosomas para devolver 25 msVCF, uno por cada cromosoma

mkdir -p "data/cromosomas_separados"
CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$data/chromosomes.txt")
if [ -z "$CHR" ]; then
  echo "No chromosome assigned to task ${SLURM_ARRAY_TASK_ID}"
  exit 0
fi

echo "Processing $CHR"
bcftools view -r "$CHR" -Oz -o "data/cromosomas_separados/${CHR}.vcf.gz" "data/msVCF_merged.vcf.gz"
bcftools index -f -t "data/cromosomas_separados/${CHR}.vcf.gz"
