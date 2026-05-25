import pandas as pd
import pysam
import re
import sys
import os

#Script para comparar dos anotaciones realizadas con distintas versiones de las dbs

arg1 = sys.argv[1] #vcf old
arg2 = sys.argv[2] #vcf new
arg3= sys.argv[3] #output


def get_clnsig(record):
    """
    Extrae el campo clínico CLNSIG desde el VCF.

    CLNSIG puede venir como:
    - string único
    - lista de valores (casos múltiples anotaciones)

    Se normaliza todo a string plano separado por '|'.
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
    Extrae genes anotados por SnpEff desde INFO/ANN.

    ANN contiene anotaciones tipo:
        allele|annotation|impact|GENE|...

    Se extrae el campo GEN (posición 4) y se eliminan duplicados.
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

def normalize_clnsig(x):
    """
    Convierte CLNSIG en un conjunto de categorías estándar.

    OBJETIVO:
    Evitar comparaciones de strings y trabajar con semántica clínica.

    Ejemplo:
        'Likely_pathogenic|benign' -> {'likely pathogenic', 'benign'}
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
    """Detecta si hay evidencia de patogenicidad."""
    return any(x in {"pathogenic", "likely pathogenic"} for x in clnsig_set)

def is_only_conflicting(clnsig_set):
    """
    Filtra variantes donde solo hay evidencia conflictiva,
    esto no se considera información clínica útil.
    """
    return clnsig_set == {"conflicting"}

def is_classified(x):
    """
    Determina si un campo CLNSIG está realmente anotado.
    Se excluyen valores vacíos o placeholders típicos de VCF.
    """
    if pd.isna(x):
        return False
    x = str(x).strip()
    return x not in ["", ".", "None", "nan"]

def vcf_to_tsv_full(vcf_file, output_tsv):
    """
    Convierte un VCF a TSV completo para análisis posterior.
    INCLUYE:
    - Campos básicos (CHROM, POS, REF, ALT)
    - Anotaciones clínicas (CLNSIG)
    - Genes (SnpEff ANN)
    - Todos los campos INFO del VCF
    - Genotipos por muestra
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
        
        # Expansión completa de campos INFO (puede ser muy grande)
        for key, val in rec.info.items():
            if isinstance(val, (tuple, list)):
                val = "|".join(map(str, val))
            row[f"INFO_{key}"] = val

        # Extracción de genotipos por muestra
        for sample in rec.samples:
            sample_data = rec.samples[sample]
            for key, val in sample_data.items():
                if isinstance(val, (tuple, list)):
                    val = ",".join(map(str, val)) 
                    
                row[f"{sample}_{key}"] = val
                
        rows.append(row)
        
    df = pd.DataFrame(rows)

    # Convertir POS a numérico para hacer el merge
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
    # Compresión gzip para eficiencia de almacenamiento
    df.to_csv(output_tsv, sep="\t", index=False, compression="gzip")
    return df

def comparar_annotaciones(old_file, new_file, output_tsv, output_reduced_tsv):
    """
    Compara dos TSV de variantes anotadas.

    OBJETIVO: Detectar cambios clínicos entre dos versiones de anotación.

    ANALIZA:
    - Variantes nuevas (solo NEW)
    - Variantes perdidas (solo OLD)
    - Variantes comunes
    - Cambios en CLNSIG
    - Ganancia/reclasificación a patogenicidad
    """
    
    old = pd.read_csv(old_file, sep="\t", dtype=str)
    new = pd.read_csv(new_file, sep="\t", dtype=str)
    
    old.columns = old.columns.str.strip()
    new.columns = new.columns.str.strip()
    
    # Clave única de variante
    key = ["CHROM", "POS", "REF", "ALT", "ID"]
    
    old["POS"] = pd.to_numeric(old["POS"], errors="coerce")
    new["POS"] = pd.to_numeric(new["POS"], errors="coerce")
    
    # Evitar duplicados antes del merge
    old = old.drop_duplicates(subset=key)
    new = new.drop_duplicates(subset=key)

    #Estadísticas de calidad de anotación 
    old_classified = old["CLNSIG"].apply(is_classified).sum()
    new_classified = new["CLNSIG"].apply(is_classified).sum()
    old_total = len(old)
    new_total = len(new)

    #Merge principal
    merged = old.merge(
        new, on=key, how="outer",
        indicator=True, suffixes=("_old", "_new")
    )
    new_variants = merged[merged["_merge"] == "right_only"]
    lost_variants = merged[merged["_merge"] == "left_only"]
    common = merged[merged["_merge"] == "both"]

    #Normalización para comparación semántica
    old_empty = common["CLNSIG_old"].apply(lambda x: not is_classified(x))
    
    old_norm = common["CLNSIG_old"].apply(normalize_clnsig)
    new_norm = common["CLNSIG_new"].apply(normalize_clnsig)
    
    changed = old_norm != new_norm

    #Detección de cambios clínicos
    now_pathogenic = new_norm.apply(is_pathogenic)
    before_pathogenic = old_norm.apply(is_pathogenic)
    
    gained_pathogenic = (~before_pathogenic) & (now_pathogenic)
    
    # Nuevas anotaciones vs reclasificación
    newly_annotated_pathogenic = old_empty & now_pathogenic
    reclassified_pathogenic = (~old_empty) & gained_pathogenic
    
    not_only_conflicting = ~new_norm.apply(is_only_conflicting)
    
    newly_annotated_pathogenic = newly_annotated_pathogenic & not_only_conflicting
    reclassified_pathogenic = reclassified_pathogenic & not_only_conflicting
    
    # Variantes de interés clínico final
    modified = common[(newly_annotated_pathogenic | reclassified_pathogenic)].copy()

    ##Seleccionar columnas dinámicamente
    # Columnas base 
    base_cols = [
        "CHROM", "POS", "REF", "ALT",
        "CLNSIG_old", "CLNSIG_new", "GENE_new"]
    # columnas de muestras (solo NEW)
    sample_cols = [col for col in modified.columns
        if col.endswith("_new") and not col.startswith("CLNSIG")]
    # Genotipos por muestra
    sample_gt_cols = [col for col in modified.columns
    if col.endswith("_GT_new")]

    #Exportar a TSV
    final_cols = [c for c in base_cols + sample_cols if c in modified.columns]
    modified = modified[final_cols].fillna(".")
    modified.to_csv(output_tsv, sep="\t", index=False)

    #Versión más reducida para visualizar mejor los datos
    final_cols_reduced = [c for c in base_cols + sample_gt_cols if c in modified.columns]
    df_reduced = modified[final_cols_reduced].fillna(".")
    df_reduced.to_csv(output_reduced_tsv, sep="\t", index=False)
    
    # LOG resumen
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

vcf_to_tsv_full(arg1, f"{arg4}/old.tsv.gz")
vcf_to_tsv_full(arg2, f"{arg4}/new.tsv.gz")
comparar_annotaciones(f"{arg4}/old.tsv.gz", f"{arg4}/new.tsv.gz",f"{arg4}/comparacion.tsv", f"{arg4}/comparacion_reducido.tsv" )
os.remove(f"{arg4}/new.tsv.gz")
os.remove(f"{arg4}/old.tsv.gz")


