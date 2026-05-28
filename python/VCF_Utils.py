import os
import subprocess
import numpy as np
import csv
import gzip
from cyvcf2 import VCF, Writer

#Este script contiene funciones para pre-procesar ficheros VCF, convirtiéndolos a matrices,comprimiéndolas y generando ficheros CSV
#temporales para la visualización de la transformación. También contiene funciones para volver a convertir estas matrices a ficheros VCF
#simplificados, separarlos en SVs y SNVs y prepararlos para la anotación con SnepEff

def vcf_to_matrix(vcf_file):
    """
    Generar matriz de genotipos a partir de un VCF multi-paciente,
    soportando variantes multialélicas y distinguiendo homocigotos y heterocigotos.
    
    Filas: variante por cada ALT, en formato CHROM:POS_ALTallelo
    Columnas: muestras
    Valores: 0 homocigoto referencia, 1 heterocigoto para ese ALT, 2 homocigoto alternativo, 3 genotipo desconocido
    """
    vcf = VCF(vcf_file)
    data = {}
    row_keys = []
    ref_list = []
    alt_list = []

    for rec in vcf:
        key = f"{rec.CHROM}:{rec.POS}" #Indicador de posicioón (cromosoma:posición)

        # Crear una columna por cada ALT. Si hay múltiples ALT, se crean varias filas para esa posición, una por cada ALT.
        for alt_index, alt_allele in enumerate(rec.ALT, start=1): 
            col_name = f"{key}_ALT{alt_index}_{alt_allele}" 
            data[col_name] = []
            ref_list.append(rec.REF)
            alt_list.append(alt_allele)
            row_keys.append(key)

        # Procesar cada muestra
        for gt in rec.genotypes:
            a1, a2 = gt[:2]  # ignorar phased si existe, da problemas al desempaquetar porque no existe para todas las muestras
            for alt_index, alt_allele in enumerate(rec.ALT, start=1):
                col_name = f"{key}_ALT{alt_index}_{alt_allele}"
                if a1 == -1 or a2 == -1:
                    val = 3  # genotipo desconocido si uno de los alelos es desconocido
                elif a1 == alt_index and a2 == alt_index:
                    val = 2  # homocigoto alternativo si ambos alelos son el mismo ALT
                elif a1 == alt_index or a2 == alt_index:
                    val = 1  # heterocigoto con uno de los alelos alternativos
                else:
                    val = 0  # homocigoto referencia
                data[col_name].append(val)

    # Convertir el diccionario a matriz numpy
    col_names = list(data.keys())
    matriz_vcf = np.array([data[k] for k in col_names], dtype='int8')

    # Asegurarse de que siempre filas = variantes, columnas = muestras
    if len(vcf.samples) == 1:
        matriz_vcf = matriz_vcf.reshape(-1, 1)
    else: # transponer para que filas sean variantes y columnas muestras
        matriz_vcf = matriz_vcf.T

    #Devuelvo la matriz, la clave y el código del paciente para poder identificarlo después al unir varias matrices
    return matriz_vcf, row_keys, ref_list, alt_list, vcf.samples


def merge_vcfs_to_npz(vcf_folder, output_npz, output_csv):
    """
    Convierte múltiples archivos VCF almacenados en una carpeta 
    a matrices y las guarda en un archivo comprimido .npz 
    """

    vcf_files = []
   # Recorremos elementos de la carpeta buscando archivos vcf
    for file in os.scandir(vcf_folder):
        if file.is_file():
            if file.name.endswith(".vcf") or file.name.endswith(".vcf.gz"):
                vcf_files.append(file.path)

    data = {} # Diccionario para almacenar matrices de variantes
    key = {}  # Diccionario para guardar samples y row_keys de cada matriz
    
    for vcf_file in vcf_files:
        # Nombre para identificar cada archivo original en la matriz npz
        id_vcf = os.path.basename(vcf_file).replace(".vcf.gz","").replace(".vcf","")
        
        # Guardar cabecera
        header_lines = []
        if vcf_file.endswith(".gz"):
            import gzip
            with gzip.open(vcf_file, 'rt') as f:
                header_lines = [line.strip() for line in f if line.startswith("#")]
        else:
            with open(vcf_file, 'r') as f:
                header_lines = [line.strip() for line in f if line.startswith("#")]

        # Convertir el VCF a matriz
        matriz, row_keys, ref_list, alt_list, samples = vcf_to_matrix(vcf_file)
        
        # Guardar todo en key
        key[f"{id_vcf}_rows"] = np.array(row_keys, dtype=object)
        key[f"{id_vcf}_ref"] = np.array(ref_list, dtype=object)
        key[f"{id_vcf}_alt"] = np.array(alt_list, dtype=object)
        key[f"{id_vcf}_samples"] = np.array(samples, dtype=object)
        key[f"{id_vcf}_header"] = np.array(header_lines, dtype=object)  #guardar cabecera
        data[id_vcf] = matriz

    # Guardar todo en .npz
    np.savez_compressed(output_npz, **data, **key)

    # Convertir a CSV global para visualización
    npz_to_csv(output_npz, output_csv)
    print(f"Guardado {len(vcf_files)} matrices en {output_npz}")


def npz_to_csv(npz_file, output_csv):
    """
    Convierte un .npz generado por merge_vcfs_to_npz
    en un único CSV global con todas las variantes y todos los pacientes.
    Rellena con 3 (genotipo desconocido) si falta información para algún paciente.
    """
    
    data = np.load(npz_file, allow_pickle=True)

    # Diccionario global: key = variante, value = dict {paciente: genotipo}
    combined = {}
    all_patients = []

    # Identificar matrices (ignorar claves auxiliares)
    matrices = [
        k for k in data.files
        if not (k.endswith("_rows") or k.endswith("_samples") or 
                k.endswith("_ref") or k.endswith("_alt") or k.endswith("_header"))
    ]

    for id_vcf in matrices:
        # Recuperar metadatos del VCF
        row_keys = data[f"{id_vcf}_rows"]
        ref_list = data[f"{id_vcf}_ref"]
        alt_list = data[f"{id_vcf}_alt"]
        samples = list(data[f"{id_vcf}_samples"])
        matriz = data[id_vcf]

        all_patients.extend(samples)

        # Construir diccionario global de variantes, pacientes y genotipos
        for i in range(len(row_keys)):
            var = f"{row_keys[i]}_{ref_list[i]}>{alt_list[i]}"
             # Identificador único para la variante (CHROM:POS_REF>ALT)
            if var not in combined:
                combined[var] = {}
            # Para cada paciente, almacenar su genotipo para esta variante
            for j, patient in enumerate(samples):
                combined[var][patient] = matriz[i, j]

    # Lista única y ordenada de pacientes
    all_patients_unique = sorted(set(all_patients))

    # Escribir CSV
    with open(output_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Variante (CHROM:POS_REF>ALT)"] + all_patients_unique)

        for var, geno_dict in combined.items():
            row = [var] + [geno_dict.get(p, 3) for p in all_patients_unique]  # rellena 3 si no hay info
            writer.writerow(row)

    print(f"CSV global exportado: {output_csv}")


def npz_to_vcf(npz_file, output_vcf):
    """
    Convierte un .npz a un VCF minimalista optimizado para anotación.
    - Una línea por variante (REF>ALT)
    - Pacientes en la columna 'PATIENTS'
    """

    data = np.load(npz_file, allow_pickle=True)

    # Identificar matrices (VCFs originales)
    matrices = [
        k for k in data.files
        if not (k.endswith("_rows") or k.endswith("_samples") or k.endswith("_ref") or k.endswith("_alt") or k.endswith("_header"))
    ]

    # Tomar cabecera del primer VCF disponible
    header_lines = []
    for id_vcf in matrices:
        header_key = f"{id_vcf}_header"
        if header_key in data.files:
            header_lines = list(data[header_key])
            break

    if not header_lines:
        # Cabecera mínima si no hay ninguna
        header_lines = [
            "##fileformat=VCFv4.2",
            '##INFO=<ID=PATIENTS,Number=.,Type=String,Description="Pacientes con la variante">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPATIENTS"
        ]
    else:
        # asegurarse de que la última línea tenga columna PATIENTS
        last_line = header_lines[-1].strip()
        cols = last_line.split("\t")
        if "PATIENTS" not in cols:
            # Si no existe, renombrar la última columna a PATIENTS
            cols[-1] = "PATIENTS"
            header_lines[-1] = "\t".join(cols)

    # Diccionario para rastrear qué pacientes portan cada variante
    variant_to_patients = {}

    # Construir diccionario 
    for id_vcf in matrices:
        matriz = data[id_vcf]
        rows = data[f"{id_vcf}_rows"]
        ref_list = data[f"{id_vcf}_ref"]
        alt_list = data[f"{id_vcf}_alt"]
        samples = data[f"{id_vcf}_samples"]

        # En cada VCF original hay solo un paciente (una muestra)
        patient = samples[0]  

        for i in range(len(rows)):
            # solo variantes presentes
            if matriz[i, 0] in [1, 2]:  
                chrom, pos = rows[i].split(":")
                ref = ref_list[i]
                alt = alt_list[i]
                var = (chrom, pos, ref, alt)
                # Registrar que este paciente porta esta variante
                if var not in variant_to_patients:
                    variant_to_patients[var] = set()
                variant_to_patients[var].add(patient)

    # Escribir VCF
    with open(output_vcf, "w") as f:
        # Escribir cabecera
        for line in header_lines:
            f.write(line + "\n")

        # Escribir variantes ordenadas
        for (chrom, pos, ref, alt), patients in sorted(variant_to_patients.items()):
            # Construir lista de pacientes separados por comas
            patient_str = ",".join(sorted(patients))
            # Columna PATIENTS en la última columna
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\t.\t{patient_str}\n")

    print(f"VCF generado: {output_vcf}")


def split_vcf(input_vcf, output_snv, output_sv):
    """
    Mueve la columna PATIENTS a INFO y elimina la columna extra final,
    añade anotaciones y divide SNV/SV.
    """
    header_lines = []
    records = []

    # Leer VCF completo
    with open(input_vcf) as f:
        for line in f:
            if line.startswith("##"):
                # Línea de metadatos de cabecera
                header_lines.append(line.strip())
            elif line.startswith("#CHROM"):
                # Línea con nombres de columnas
                header_lines.append(line.strip())
                columns = line.strip().split("\t")
                # Encontrar índice de la columna PATIENTS
                patients_col_idx = columns.index("PATIENTS")

            # Ajustar la línea #CHROM para eliminar PATIENTS
                header_lines_cleaned = []
                for h in header_lines:
                    if h.startswith("#CHROM"):
                        cols = h.split("\t")
                        # Eliminar columna PATIENTS de la cabecera
                        cols = cols[:patients_col_idx] + cols[patients_col_idx+1:]
                        header_lines_cleaned.append("\t".join(cols))
                    else:
                        header_lines_cleaned.append(h)            
            else:
                records.append(line.strip().split("\t"))

    snv_lines = []
    sv_lines = []

    # Añadir INFO header para PATIENTS y ANNOT
    info_header_patients = '##INFO=<ID=PATIENTS,Number=.,Type=String,Description="Pacientes con la variante">'
    info_header_annot = '##INFO=<ID=ANNOT,Number=1,Type=String,Description="Anotación extra">'
    header_lines.insert(-1, info_header_patients)
    header_lines.insert(-1, info_header_annot)

    # Procesar registros
    for rec in records:
        ref = rec[3]
        alt = rec[4]
        pacientes = rec[patients_col_idx]

        # Construir INFO
        info_field = f'PATIENTS={pacientes};'
        rec[7] = info_field  # INFO es columna 7

        # Eliminar columna PATIENTS original
        rec_cleaned = rec[:patients_col_idx] + rec[patients_col_idx+1:]

        # Separar SNV vs SV
        if len(ref) == 1 and all(len(a) == 1 for a in alt.split(",")):
            snv_lines.append("\t".join(rec_cleaned))
        else:
            sv_lines.append("\t".join(rec_cleaned))

    # Escribir VCFs de salida 
    with open(output_snv, "w") as f:
        for h in header_lines:
            f.write(h + "\n")
        for r in snv_lines:
            f.write(r + "\n")

    with open(output_sv, "w") as f:
        for h in header_lines:
            f.write(h + "\n")
        for r in sv_lines:
            f.write(r + "\n")


def compress_vcf(input_vcf, output_vcf_gz):
    """
    Comprime un VCF en formato .gz usando gzip de Python.
    """

    # Leer el archivo en chunks para optimizar memoria con archivos grandes
    with open(input_vcf, 'rb') as f_in, gzip.open(output_vcf_gz, 'wb') as f_out:
        # Procesar en bloques de 1MB
        for chunk in iter(lambda: f_in.read(1024*1024), b""):
            f_out.write(chunk)
    print(f"Archivo comprimido: {output_vcf_gz}")

def prepare_vcf(input_vcf, provisional_vcf):
    """
    Prepara un VCF para SnpEff:
    - Mueve la columna PATIENTS a INFO como MY_PATIENTS
    - No añade MY_PATIENTS si no hay valor
    - Elimina la columna PATIENTS de la tabla
    """
    header_lines = []
    records = []

    # Procesar línea por línea
    with open(input_vcf) as f:
        for line in f:
            if line.startswith("##"):
                # Ajusta cabecera INFO si existe PATIENTS
                if line.startswith('##INFO=<ID=PATIENTS'):
                    info_line = line.replace('ID=PATIENTS', 'ID=MY_PATIENTS') \
                                    .replace('Pacientes con la variante', 'Pacientes con la variante (renombrado para SnpEff)')
                    header_lines.append(info_line)
                else:
                    header_lines.append(line.strip())
            elif line.startswith("#CHROM"):
                # Línea con nombres de columnas
                cols = line.strip().split("\t")
                # Encontrar el índice de la columna PATIENTS y eliminarla de la cabecera
                patients_col_idx = cols.index("PATIENTS")
                header_lines.append("\t".join(cols[:patients_col_idx] + cols[patients_col_idx+1:]))
            else:
                fields = line.strip().split("\t")
                # Obtener valor de la columna PATIENTS
                if len(fields) <= patients_col_idx:
                    pacientes_value = ""
                else:
                    pacientes_value = fields[patients_col_idx]

                # Construir INFO
                info = fields[7] if fields[7] != "." else ""
                if pacientes_value not in ("", "."):
                    if info:
                        info += f";MY_PATIENTS={pacientes_value}"
                    else:
                        info = f"MY_PATIENTS={pacientes_value}"
                fields[7] = info

                # Eliminar columna PATIENTS
                if len(fields) > patients_col_idx:
                    fields = fields[:patients_col_idx] + fields[patients_col_idx+1:]

                records.append("\t".join(fields))

    # Escribir archivo VCF preparado
    with open(provisional_vcf, "w") as f:
        for h in header_lines:
            f.write(h + "\n")
        for r in records:
            f.write(r + "\n")



