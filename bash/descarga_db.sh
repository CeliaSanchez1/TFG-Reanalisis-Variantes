#!/bin/bash
#SBATCH --job-name=descargar_db    # Job name
#SBATCH --output=descargar_db_%j.out       # Output file (%j expands to job ID)
#SBATCH --error=descargar_db_%j.err        # Error file
#SBATCH --time=02:00:00             # Walltime (3 minute)
#SBATCH --ntasks=1                  # Number of tasks (processes)
#SBATCH --cpus-per-task=4           # Number of CPU cores per task
#SBATCH --mem=8G                 # Memory per node
#SBATCH --partition=         # Queue/partition (adjust to your system)

#Script para comprobar las versiones de las db y volver a descargarlas si no son las más recientes 
ml Java

if [ ! -f "data/snpEff/snpEff.config" ]; then
  echo "data/snepEff/snpEff.config not found"; exit
fi

WORKDIR="snpEff/dbs" 
mkdir -p "$WORKDIR"

#Descarga de clinVar
ClinVarURL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/" #URL REPOSITORIO DB
clinvar_version_file="${WORKDIR}/clinvar.version" #archivo de version actual
clinvarFile=$(curl -s "$ClinVarURL/" | grep -oE 'clinvar_[0-9]+\.vcf\.gz' | sort -V | tail -n1) #Busqueda de la version en el repositorio
clinvar_path="${WORKDIR}/${clinvarFile}" 
clinvar_tbi="${clinvar_path}.tbi"

#Comprobar si la version del repositorio coincide con la descargada y si no descargarla
if [[ -f "$clinvar_version_file" ]] && grep -qx "$clinvarFile" "$clinvar_version_file" \
   && [[ -s "$clinvar_path" ]] && [[ -s "$clinvar_tbi" ]]; then
  echo "ClinVar up to date" 
else 
  echo "Downloading ClinVar $clinvarFile"
  cd "$WORKDIR"
  rm -f clinvar_*.vcf.gz clinvar_*.vcf.gz.tbi
  wget -q -O "$clinvarFile" "$ClinVarURL/$clinvarFile"
  wget -q -O "${clinvarFile}.tbi" "$ClinVarURL/${clinvarFile}.tbi"
  if [[ ! -s "$clinvarFile" ]] || [[ ! -s "${clinvarFile}.tbi" ]]; then
    echo "ERROR: ClinVar download incomplete or corrupted"
    exit 1
  fi
  clinvar_path="${WORKDIR}/${clinvarFile}"
  echo "$clinvarFile" > "$clinvar_version_file"
  REANNOTATE=true
  echo "ClinVar ready"
fi

# dbSNP
dbsnpURL="https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/"
dbsnp_version_file="${WORKDIR}/dbsnp.version"
dbsnpFile=$(curl -s "$dbsnpURL/" | grep -oE '[0-9]+-All\.vcf\.gz' | sort -V | tail -n1)
if [[ -z "$dbsnpFile" ]]; then
  echo "ERROR: dbSNP scraping failed (no file detected)"
  exit 1
fi
dbsnp_path="${WORKDIR}/${dbsnpFile}"
dbsnp_tbi="${dbsnp_path}.tbi"
echo "Detected dbSNP version: $dbsnpFile"
if [[ -f "$dbsnp_version_file" ]] && grep -qx "$dbsnpFile" "$dbsnp_version_file"; then
  echo "dbSNP up to date"
else
  echo "New dbSNP version detected"
  cd "$WORKDIR"
  rm -f GCF_*.gz GCF_*.gz.tbi
  wget -q -O "$dbsnpFile" "$dbsnpURL/$dbsnpFile"
  wget -q -O "${dbsnpFile}.tbi" "$dbsnpURL/${dbsnpFile}.tbi"
  if [[ ! -s "$dbsnpFile" ]]; then
    echo "ERROR: dbSNP VCF download failed or empty"
    exit 1
  fi
  if [[ ! -s "${dbsnpFile}.tbi" ]]; then
    echo "ERROR: dbSNP index download failed"
    exit 1
  fi
  echo "$dbsnpFile" > "$dbsnp_version_file"
  REANNOTATE=true
  echo "dbSNP updated successfully"
fi

# dbNSFP
dbnsfpURL="http://dbnsfp.houstonbioinformatics.org/dbNSFPzip"
dbnsfp_version_file="${WORKDIR}/dbnsfp.version"
mkdir -p "$WORKDIR"
dbnsfpFile=$(curl -s "$dbsnfpURL/" | grep -oE 'dbNSFP4+-[0-9]\.txt\.gz' | sort -V | tail -n1)
dbnsfp_path="${WORKDIR}/${dbnsfpFile}"
echo "Using dbNSFP file: $dbnsfpFile"
if [[ -f "$dbnsfp_version_file" ]] && grep -qx "$dbnsfpFile" "$dbnsfp_version_file"; then
  echo "dbNSFP up to date"
  exit 0
fi
echo "Downloading dbNSFP..."
cd "$WORKDIR"
rm -f "$dbnsfpFile"
wget -O "$dbnsfpFile"  "$dbnsfpURL/$dbnsfpFile"
if [[ ! -s "$dbnsfpFile" ]]; then
  echo "ERROR: dbNSFP download failed (empty or blocked)"
  exit 1
fi
echo "$dbnsfpFile" > "$dbnsfp_version_file"
REANNOTATE=true
echo "dbNSFP updated successfully"

#Descargar MITOMAP para anotar el cromosoma M
wget -O data/mitomap_disease.vcf "https://mitomap.org/cgi-bin/disease.cgi?format=vcf"
