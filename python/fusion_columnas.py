"""
Fusiona varios TSV (mismas columnas, filas distintas) en uno solo.
La cabecera se escribe una sola vez desde el primer archivo.

Uso:
    python merge_tsvs.py -i *.changed.tsv -o merged.tsv
    python merge_tsvs.py -i chr1.tsv chr2.tsv chr3.tsv -o merged.tsv
"""
import argparse
import glob
import os
import sys

def merge(inputs: list, output: str):
    # Expandir globs por si se pasa "*.tsv" como string literal
    files = []
    for pattern in inputs:
        expanded = glob.glob(pattern)
        files.extend(sorted(expanded) if expanded else [pattern])

    if not files:
        sys.exit("Error: no se encontraron archivos de entrada.")
    header = None
    total_rows = 0
    with open(output, "w") as out:
        for path in files:
            if not os.path.isfile(path):
                print(f"  [WARN] No existe: {path}", file=sys.stderr)
                continue
            with open(path) as f:
                file_header = f.readline()
                if header is None:
                    header = file_header
                    out.write(header)
                elif file_header != header:
                    sys.exit(
                        f"Error: las columnas de '{path}' no coinciden con las del primer archivo.\n"
                        f"  Esperado: {header.strip()}\n"
                        f"  Obtenido: {file_header.strip()}"
                    )

                rows = 0
                for line in f:
                    out.write(line)
                    rows += 1

            print(f"  {path}: {rows} variantes")
            total_rows += rows

    print(f"\nTotal: {total_rows} variantes → {output}")


def main():
    p = argparse.ArgumentParser(description="Fusiona TSVs con las mismas columnas")
    p.add_argument("-i", "--inputs", nargs="+", required=True,
                   metavar="FILE", help="Archivos TSV de entrada (acepta globs)")
    p.add_argument("-o", "--output", required=True,
                   metavar="FILE", help="Archivo TSV de salida")
    args = p.parse_args()
    merge(args.inputs, args.output)

if __name__ == "__main__":
    main()