import pandas as pd
import sys

arg1 = sys.argv[7] #RUTA TABLA RELACIÓN HPO-GENES FORMATO TXT
arg2 = sys.argv[8] #Carpeta para el tsv de salida para cada gen
arg3 = sys.argv[9] #Carpeta para el tsv de salida reducido para cada gen
arg4= sys.argv[10]#carpeta tsv final salida
arg5= sys.argv[11]#carpeta tsv final salida reducido
arg6= sys.argv[12] #Directorio de salida

##Introduzco la información de los HPO a la tabla de variantes modificadas
var = pd.read_csv(arg8 sep="\t")
var_reduced = pd.read_csv(arg9 , sep="\t")
hpo = pd.read_csv(arg7, sep="\t")

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
merged.to_csv(arg10, sep="\t", index=False)
merged_reduced.to_csv(arg11, sep="\t", index=False)

os.remove(arg9)
os.remove(arg8)