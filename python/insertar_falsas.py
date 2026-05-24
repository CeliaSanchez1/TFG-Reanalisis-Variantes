"""
insert_fake_variants.py
-----------------------
Inserta variantes falsas en una copia del VCF de MITOMAP para testear
la función de comparación de versiones de anotación.
"""

import argparse
import os
import sys

# Variantes falsas a insertar
FAKE_VARIANTS = [
    ("MT", 10,  "FAKE_INSERTED_FOR_COMPARISON_TESTING", "T",  "C",      ".", ".", "AC=2;AF=0.001;aachange=noncoding;homoplasmy=nr;heteroplasmy=nr;Disease=FAKE;DiseaseStatus=reported-lp"),
    ("MT", 56,  "FAKE_INSERTED_FOR_COMPARISON_TESTING", "A",  "AC,ATC", ".", ".", "AC=1;AF=0.001;aachange=noncoding;homoplasmy=nr;heteroplasmy=nr;Disease=FAKE;DiseaseStatus=cfrm-p"),
]

def parse_args():
    parser = argparse.ArgumentParser(
        description="Inserta variantes falsas en el VCF de MITOMAP."
    )
    parser.add_argument("--input",  "-i", required=True, help="VCF original de MITOMAP.")
    parser.add_argument("--output", "-o", required=True, help="VCF de salida con variantes falsas.")
    return parser.parse_args()

def main():
    args = parse_args()

    if not os.path.isfile(args.input):
        print(f"[ERROR] No se encuentra el archivo: {args.input}", file=sys.stderr)
        sys.exit(1)

    header_lines = []
    data_lines   = []

    with open(args.input) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                header_lines.append(line)
            elif line.strip():
                data_lines.append(line)

    fake_lines = [
        "\t".join([chrom, str(pos), id_, ref, alt, qual, filter_, info])
        for chrom, pos, id_, ref, alt, qual, filter_, info in FAKE_VARIANTS
    ]

    all_lines = data_lines + fake_lines
    all_lines.sort(key=lambda l: int(l.split("\t")[1]))

    with open(args.output, "w") as fh:
        for h in header_lines:
            fh.write(h + "\n")
        for line in all_lines:
            fh.write(line + "\n")

    print(f"[INFO] Variantes originales : {len(data_lines)}")
    print(f"[INFO] Variantes falsas      : {len(fake_lines)}")
    print(f"[INFO] Total en salida       : {len(all_lines)}")
    print(f"[INFO] Archivo escrito en    : {args.output}")

if __name__ == "__main__":
    main()