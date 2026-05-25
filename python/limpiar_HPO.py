import pandas as pd
import sys

#Script para limpiar la tabla introducir la información de los HPO al TSV de variantes

arg1 = sys.argv[1] #RUTA TABLA RELACIÓN HPO-GENES FORMATO TXT
arg2 = sys.argv[2] #Carpeta para el tsv de salida para cada gen
arg3 = sys.argv[3] #Carpeta para el tsv de salida reducido para cada gen
arg4= sys.argv[4]#carpeta tsv final salida
arg5= sys.argv[5]#carpeta tsv final salida reducido
arg6= sys.argv[6] #Directorio de salida

##Introduzco la información de los HPO a la tabla de variantes modificadas
var = pd.read_csv(arg2 sep="\t")
var_reduced = pd.read_csv(arg3 , sep="\t")
hpo = pd.read_csv(arg1, sep="\t")

#Elimino filas duplicadas basándome en las columnas "gene_symbol" y "hpo_id", para asegurar que cada combinación de gen y término HPO sea única
hpo = hpo.drop_duplicates(subset=["gene_symbol", "hpo_id"])
#Agrupo los términos HPO por gen y los uno en una sola cadena separada por "|", para evitar tener múltiples entradas por variante
hpo_agg = hpo.groupby("gene_symbol")["hpo_id"].apply(lambda x: "|".join(sorted(set(x)))).reset_index()

#Divido la columna "GENE_new" en varias filas para manejar casos donde hay múltiples genes asociados a una variante
var["GENE_new"] = var["GENE_new"].str.split("|")
var = var.explode("GENE_new")
var_reduced["GENE_new"] = var_reduced["GENE_new"].str.split("|")
var_reduced = var_reduced.explode("GENE_new")

merged = var.merge(hpo_agg, left_on="GENE_new", right_on="gene_symbol", how="left").drop(columns=["gene_symbol"])
merged_reduced = var_reduced.merge(hpo_agg, left_on="GENE_new", right_on="gene_symbol", how="left").drop(columns=["gene_symbol"])
merged.to_csv(arg4, sep="\t", index=False)
merged_reduced.to_csv(arg5, sep="\t", index=False)

os.remove(arg3)
os.remove(arg2)
