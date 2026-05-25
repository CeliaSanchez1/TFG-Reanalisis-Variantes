import argparse
import os
import sys

#Script para insertar un par de variantes falsas a la base de datos de MITOMAP para comprobar la detección de las variantes nuevas

# Variantes a insertar
FAKE_VARIANTS = [
    ("MT", 10,  "FAKE_INSERTED_FOR_COMPARISON_TESTING", "T",  "C",      ".", ".", "AC=2;AF=0.001;aachange=noncoding;homoplasmy=nr;heteroplasmy=nr;Disease=FAKE;DiseaseStatus=reported-lp"),
    ("MT", 56,  "FAKE_INSERTED_FOR_COMPARISON_TESTING", "A",  "AC,ATC", ".", ".", "AC=1;AF=0.001;aachange=noncoding;homoplasmy=nr;heteroplasmy=nr;Disease=FAKE;DiseaseStatus=cfrm-p"),
]

def parse_args():
    """
    Definir y procesar argumentos de línea de comandos.
    """

    # Crea el parser principal del programa
    parser = argparse.ArgumentParser(
        description="Inserta variantes falsas en el VCF de MITOMAP."
    )
    # Archivo VCF de entrada
    parser.add_argument("--input",  "-i", required=True, help="VCF original de MITOMAP.")
    # Archivo VCF de salida
    parser.add_argument("--output", "-o", required=True, help="VCF de salida con variantes falsas.")
    return parser.parse_args()

def main():
    """
    Función para leer argumentos CLI, separar cabeceras, generar variantes falsas y añadirlas y reordenar el VCF. 
    """

    # Obtiene argumentos de línea de comandos
    args = parse_args()
    if not os.path.isfile(args.input):
        print(f"[ERROR] No se encuentra el archivo: {args.input}", file=sys.stderr)
        sys.exit(1)

    header_lines = [] # Lista para cabecera
    data_lines   = [] # Lista para variantes

    with open(args.input) as fh: 
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                #Extraer líneas de cabecera
                header_lines.append(line)
            elif line.strip():
                #Extraer líneas de datos
                data_lines.append(line)

    #Prepara las variantes falsas en formato VCF
    fake_lines = [
        "\t".join([chrom, str(pos), id_, ref, alt, qual, filter_, info])
        for chrom, pos, id_, ref, alt, qual, filter_, info in FAKE_VARIANTS
    ]
    
    #Une variantes propias de la db con las falsas y las ordena por posición genómica
    all_lines = data_lines + fake_lines
    all_lines.sort(key=lambda l: int(l.split("\t")[1]))

    #Escribe el nuevo VCF
    with open(args.output, "w") as fh:
        for h in header_lines:
            #Primero cabecera
            fh.write(h + "\n")
        for line in all_lines:
            #Luego datos
            fh.write(line + "\n")

    print(f"[INFO] Variantes originales : {len(data_lines)}")
    print(f"[INFO] Variantes falsas      : {len(fake_lines)}")
    print(f"[INFO] Total en salida       : {len(all_lines)}")
    print(f"[INFO] Archivo escrito en    : {args.output}")

if __name__ == "__main__":
    main()
