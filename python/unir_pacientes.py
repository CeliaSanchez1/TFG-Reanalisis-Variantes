import pandas as pd

#Script para unir la información clínica de los pacientes en el archivo de variantes a través de los HPO

variants = pd.read_csv("/data/merged.changed.tsv", sep="\t", header=0)
patients  = pd.read_csv("/data/HPO_Pacientes.tsv", sep="\t", header=0)
fixed_cols = ["VARIANT", "OLD", "NEW", "HPO"]
sample_cols = [c for c in variants.columns if c not in fixed_cols]

#Reconocer solo genotipos válidos y distintos de homocigoto de referencia
long = (
    variants
    .melt(id_vars=fixed_cols, value_vars=sample_cols,
          var_name="Local_ID", value_name="genotype")
    .query("genotype in ['0/1', '0/2', '1/0','1/1', '1/2', '2/0', '2/1', '2/2']")
)

#Join con tabla de pacientes por Local_ID
merged = long.merge(patients[["Local_ID", "HPO_Observed_IDs"]], on="Local_ID")

#Calcular intersección y guardar solo los HPO coincidentes
def hpo_match(row):
    gene_hpos    = set(row["HPO"].split("|"))
    patient_hpos = set(row["HPO_Observed_IDs"].split(";"))
    shared       = gene_hpos & patient_hpos
    return ";".join(shared) if shared else None

merged["HPO_match"] = merged.apply(hpo_match, axis=1)

#Filtrar solo variantes relevantes y columnas de interés 
result = (
    merged[merged["HPO_match"].notna()]
    [["VARIANT", "OLD", "NEW", "Local_ID", "genotype", "HPO_match"]]
    .copy()
)

# Extraer cromosoma y posición para ordenar
result["chrom"] = result["VARIANT"].str.split(":").str[0].str.replace("chr", "")
result["pos"]   = result["VARIANT"].str.split(":").str[1].astype(int)

# Orden cromosómico correcto (1-22, X, Y, MT)
chrom_order = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
result["chrom"] = pd.Categorical(result["chrom"], categories=chrom_order, ordered=True)
result = result.sort_values(["chrom", "pos"]).drop(columns=["chrom", "pos"])

#Exportar
result.to_csv("/data/variantes_relevantes.tsv", sep="\t", index=False)
