import argparse
import glob
import os
import sys

#Script para fusionar varios TSV (mismas columnas, filas distintas) en uno solo
#Uso python merge_tsvs.py -i *.changed.tsv -o merged.tsv

def merge(inputs: list, output: str):
    """
    Fusiona múltiples archivos TSV en un único archivo de salida
    """

    # Expandir globs por si se pasa "*.tsv" como string literal
    files = []
    for pattern in inputs:
        # Busca coincidencias del patrón
        expanded = glob.glob(pattern)
        # Si encuentra archivos, los añade ordenados
        files.extend(sorted(expanded) if expanded else [pattern])

    if not files:
        sys.exit("Error: no se encontraron archivos de entrada.")
        
    header = None
    total_rows = 0
    
    with open(output, "w") as out:
        #Recorre los archivos TSV
        for path in files:
            if not os.path.isfile(path):
                print(f"  [WARN] No existe: {path}", file=sys.stderr)
                continue
            with open(path) as f:
                #Recorre el archivo
                # Lee la primera línea (cabecera)
                file_header = f.readline()
                #Guarda la cabecera solo del primer TSV 
                if header is None:
                    header = file_header
                    #Escribe la cabecera en el archivo final
                    out.write(header)
                elif file_header != header:
                    #Si alguno de los TSV contiene una cabecera diferente devuelve una advertencia
                    sys.exit(
                        f"Error: las columnas de '{path}' no coinciden con las del primer archivo.\n"
                        f"  Esperado: {header.strip()}\n"
                        f"  Obtenido: {file_header.strip()}"
                    )

                rows = 0
                for line in f:
                    #Copia cada línea del archivo detrás de la cabecera
                    out.write(line)
                    rows += 1

            print(f"  {path}: {rows} variantes")
            total_rows += rows
    print(f"\nTotal: {total_rows} variantes en {output}")


def main():
    # Crea parser de argumentos CLI
    p = argparse.ArgumentParser(description="Fusiona TSVs con las mismas columnas")
     # Argumento de archivos de entrada
    p.add_argument("-i", "--inputs", nargs="+", required=True,
                   metavar="FILE", help="Archivos TSV de entrada (acepta globs)")
    # Argumento de archivo de salida
    p.add_argument("-o", "--output", required=True,
                   metavar="FILE", help="Archivo TSV de salida")
    # Parsea argumentos CLI y ejecuta la función
    args = p.parse_args()
    merge(args.inputs, args.output)

if __name__ == "__main__":
    main()
