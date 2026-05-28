import os
import subprocess
import numpy as np
import csv
import gzip
import sys
from cyvcf2 import VCF, Writer

#Script para anotar archivo VCF empleando SnpEff. Se emplea con el Wrapper correspondiente

arg1 = sys.argv[1] #INPUT A ANOTAR
arg2 = sys.argv[2] #OUTPUT ANOTADO Y COMPRIMIDO

def annotate_snv_SnpEff(input_vcf, output_vcf_gz, genoma_referencia):
    """
    Ejecuta SnpEff y comprime a .gz 
    """

    #Ejecución de snpEff como proceso externo, capturando salida, errores y mensajes de diagnóstico
    cmd = ["snpEff", "-v", genoma_referencia, input_vcf]
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    #Comprimir el resultado en formato gz
    with gzip.open(output_vcf_gz, "wt") as out_f:
        for line in proc.stdout:
            out_f.write(line)
            
    stderr = proc.stderr.read()
    return_code = proc.wait()
    
    #Captura posibles errores
    if return_code != 0:
        print("ERROR en SnpEff:")
        print(stderr)
        raise RuntimeError("SnpEff falló")
    print(f"VCF anotado comprimido generado: {output_vcf_gz}")

annotate_snv_SnpEff(arg1, arg2, "GRCh38.99")
