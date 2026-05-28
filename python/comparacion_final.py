import pandas as pd
import pysam
import re
import sys
import os

#Script para comparar anotaciones con distintas versiones de las dbs sobre un mismo VCF. Convierte VCF a TSV, extrayendo campos
#relevantes para el diagnóstico, y compara TSVs entr esí recogiendo las variantes que pasan a clasificación de patogenicidad. 
#Genera un log con las estadísticas de la comparación y un TSV con las variantes relevantes.
#Se añade a este script también la incorporación de la información de los HPO asociados

arg1 = sys.argv[1] #txt relación HPO-paciente
arg2 = sys.argv[2] #VCF antiguo
arg3= sys.argv[3] #VCF nuevo
arg6= sys.argv[4]  #Directorio de salida
arg10= sys.argv[5]#carpeta tsv final salida
arg11= sys.argv[6]#carpeta tsv final salida reducido


##FUNCIONES PARA PROCESAR VCF COMO TSV

def get_clnsig(record):
     """
    Extrae el campo CLNSIG (significado clínico) desde un registro VCF.

    Este campo puede contener múltiples valores (lista/tuple), que se
    concatenan con '|'. Si no existe, devuelve None.
    """
    try:
        val = record.info.get("CLNSIG")
        if val is None:
            return None
        if isinstance(val, (tuple, list)):
            return "|".join(map(str, val))
        return str(val)
    except:
        return None

def get_gene_snpeff(record):
    """
    Extrae los genes anotados por SnpEff desde el campo INFO/ANN.
    ANN contiene anotaciones separadas por '|'. El gen se encuentra en la
    posición 4 (índice 3). Se devuelven genes únicos separados por '|'.
    """
    try:
        ann = record.info.get("ANN")
        if ann is None:
            return None
        genes = set()
        for entry in ann:
            fields = entry.split("|")
            if len(fields) > 3:
                gene = fields[3]
                if gene:
                    genes.add(gene)
        
        return "|".join(sorted(genes)) if genes else None
    except:
        return None


##FUNCIONES PARA COMPARAR VERSIONES DE TSV

def normalize_clnsig(x):
    """
    Normaliza el campo CLNSIG para facilitar comparaciones entre versiones.
    - Convierte a minúsculas
    - Divide por '|', ',', ';'
    - Agrupa términos en categorías:
        pathogenic, likely pathogenic, benign, uncertain, conflicting
    """
    if pd.isna(x):
        return set()
    import re
    parts = re.split(r"[|,;]", str(x).lower())
    normalized = set()
    for p in parts:
        p = p.strip().replace("_", " ")   
        if "conflicting" in p:
            normalized.add("conflicting")
        elif "likely" in p and "pathogenic" in p:
            normalized.add("likely pathogenic")
        elif "pathogenic" in p:
            normalized.add("pathogenic")
        elif "benign" in p:
            normalized.add("benign")
        elif "uncertain" in p:
            normalized.add("uncertain")
        elif p:
            normalized.add(p)
    return normalized

def is_pathogenic(clnsig_set):
    """
    Devuelve True si el conjunto CLNSIG contiene evidencia patogénica.
    """
    return any(x in {"pathogenic", "likely pathogenic"} for x in clnsig_set)

def is_only_conflicting(clnsig_set):
    """
    Devuelve True si el único valor es 'conflicting'.
    """
    return clnsig_set == {"conflicting"}

def is_classified(x):
    """
    Indica si un CLNSIG está correctamente clasificado (no vacío ni nulo).
    """
    if pd.isna(x):
        return False
    x = str(x).strip()
    return x not in ["", ".", "None", "nan"]


def vcf_to_tsv_full(vcf_file, output_tsv):
    """
    Convierte un archivo VCF en un TSV enriquecido.

    Extrae:
    - Campos básicos del VCF
    - CLNSIG
    - Genes (SnpEff)
    - Todos los campos INFO
    - Genotipos de muestras

    Guarda el resultado en formato TSV comprimido (gzip).
    """
    vcf = pysam.VariantFile(vcf_file)
    rows = []
    for rec in vcf.fetch():
        row = {
            "CHROM": rec.chrom,
            "POS": rec.pos,
            "ID": rec.id,
            "REF": rec.ref,
            "ALT": ",".join(rec.alts) if rec.alts else None,
            "QUAL": rec.qual,
            "FILTER": ";".join(rec.filter.keys()) if rec.filter else None,
            "CLNSIG": get_clnsig(rec),
            "GENE": get_gene_snpeff(rec)
        }
        
        # Campos INFO del VCF
        for key, val in rec.info.items():
            if isinstance(val, (tuple, list)):
                val = "|".join(map(str, val))
            row[f"INFO_{key}"] = val

        # Campos de muestras (genotipos y otros)
        for sample in rec.samples:
            sample_data = rec.samples[sample]
            for key, val in sample_data.items():
                if isinstance(val, (tuple, list)):
                    val = ",".join(map(str, val))  
                row[f"{sample}_{key}"] = val
        rows.append(row)
        
    df = pd.DataFrame(rows)
    # Asegurar tipo numérico en POS
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
    df.to_csv(output_tsv, sep="\t", index=False, compression="gzip")
    return df


def comparar_annotaciones(old_file, new_file, output_tsv, output_reduced_tsv):
    """
    Compara dos TSV de variantes (OLD vs NEW) y detecta:

    - Variantes nuevas
    - Variantes perdidas
    - Variantes comunes
    - Cambios en CLNSIG
    - Ganancia de anotaciones patogénicas

    Genera:
    - TSV completo de cambios
    - TSV reducido (para visualización)
    - Log resumen
    """
    
    old = pd.read_csv(old_file, sep="\t", dtype=str)
    new = pd.read_csv(new_file, sep="\t", dtype=str)
    
    old.columns = old.columns.str.strip()
    new.columns = new.columns.str.strip()
    
    key = ["CHROM", "POS", "REF", "ALT", "ID"]
    
    old["POS"] = pd.to_numeric(old["POS"], errors="coerce")
    new["POS"] = pd.to_numeric(new["POS"], errors="coerce")
    
    old = old.drop_duplicates(subset=key)
    new = new.drop_duplicates(subset=key)

    # Estadísticas de clasificación
    old_classified = old["CLNSIG"].apply(is_classified).sum()
    new_classified = new["CLNSIG"].apply(is_classified).sum()
    
    old_total = len(old)
    new_total = len(new)

    # Merge para comparar variantes
    merged = old.merge(
        new, on=key, how="outer",
        indicator=True, suffixes=("_old", "_new")
    )
    
    new_variants = merged[merged["_merge"] == "right_only"]
    lost_variants = merged[merged["_merge"] == "left_only"]
    common = merged[merged["_merge"] == "both"]

    # CLNSIG processing
    old_empty = common["CLNSIG_old"].apply(lambda x: not is_classified(x))
    
    old_norm = common["CLNSIG_old"].apply(normalize_clnsig)
    new_norm = common["CLNSIG_new"].apply(normalize_clnsig)
    
    changed = old_norm != new_norm
    
    now_pathogenic = new_norm.apply(is_pathogenic)
    before_pathogenic = old_norm.apply(is_pathogenic)
    
    gained_pathogenic = (~before_pathogenic) & (now_pathogenic)
    newly_annotated_pathogenic = old_empty & now_pathogenic
    reclassified_pathogenic = (~old_empty) & gained_pathogenic

    # Filtrar solo casos relevantes
    not_only_conflicting = ~new_norm.apply(is_only_conflicting)
    
    newly_annotated_pathogenic = newly_annotated_pathogenic & not_only_conflicting
    reclassified_pathogenic = reclassified_pathogenic & not_only_conflicting
    
    modified = common[(newly_annotated_pathogenic | reclassified_pathogenic)].copy()

    # seleccionar columnas dinámicamente
    base_cols = [
        "CHROM", "POS", "REF", "ALT",
        "CLNSIG_old", "CLNSIG_new", "GENE_new"]

    # columnas de muestras (solo NEW)
    sample_cols = [col for col in modified.columns
        if col.endswith("_new") and not col.startswith("CLNSIG")
    ]

    sample_gt_cols = [col for col in modified.columns
    if col.endswith("_GT_new")]

    #Versión completa
    final_cols = base_cols + sample_cols
    final_cols = [c for c in final_cols if c in modified.columns]
    modified = modified[final_cols].fillna(".")
    modified.to_csv(output_tsv, sep="\t", index=False)
    
    #Versión más reducida para visualizar mejor los datos
    final_cols_reduced = base_cols + sample_gt_cols
    final_cols_reduced = [c for c in final_cols_reduced if c in modified.columns]
    df_reduced = modified[final_cols_reduced].fillna(".")
    df_reduced.to_csv(output_reduced_tsv, sep="\t", index=False)

    # LOG de resumen
    with open("comparison_summary.log", "w") as log:
        log.write("RESUMEN DE COMPARACIÓN DE ANOTACIONES\n\n")
        log.write(f"Variantes en OLD: {len(old)}\n")
        log.write(f"Variantes en NEW: {len(new)}\n\n")
        log.write(f"Variantes en ambas anotaciones: {len(common)}\n")
        log.write(f"Variantes nuevas (solo NEW): {len(new_variants)}\n")
        log.write(f"Variantes perdidas (solo OLD): {len(lost_variants)}\n\n")
        log.write("CLNSIG coverage:\n")
        log.write(f"OLD clasificados: {old_classified}/{old_total} "
                  f"({old_classified/old_total*100:.2f}%)\n")
        log.write(f"NEW clasificados: {new_classified}/{new_total} "
                  f"({new_classified/new_total*100:.2f}%)\n\n")
        log.write(f"Variantes nuevas patogénicas (antes sin CLNSIG): {newly_annotated_pathogenic.sum()}\n")
        log.write(f"Variantes reclasificadas a patogénicas: {reclassified_pathogenic.sum()}\n")
    return common, modified

vcf_to_tsv_full(arg2, f"{arg4}/old.tsv.gz")
vcf_to_tsv_full(arg3, f"{arg4}/new.tsv.gz")
comparar_annotaciones(f"{arg4}/old.tsv.gz", f"{arg4}/new.tsv.gz", f"{arg4}/comparacion_anotaciones.tsv", f"{arg4}/comparacion_anotaciones_reducido.tsv")
os.remove(f"{arg4}/new.tsv.gz")
os.remove(f"{arg4}/old.tsv.gz")

##Integración con información de fenotipo
var = pd.read_csv(f"{arg4}/comparacion_anotaciones.tsv", sep="\t")
var_reduced = pd.read_csv(f"{arg4}/comparacion_anotaciones_reducido.tsv", sep="\t")
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
merged.to_csv(arg5, sep="\t", index=False)
merged_reduced.to_csv(arg6, sep="\t", index=False)

os.remove(f"{arg4}/comparacion_anotaciones_reducido.tsv")
os.remove(f"{arg4}/comparacion_anotaciones.tsv" )

