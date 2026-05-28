import os
import numpy as np
from cyvcf2 import VCF #Librería para trabajar con archivos de llamado de variantes
from VCF_Utils import merge_vcfs_to_npz, npz_to_vcf, split_vcf 

#Script que contiene el flujo de trabajo de creación de matrices npz, reconversión a VCF y separación en SV y SNV

merge_vcfs_to_npz(r"/data/vcf_original", r"/data/pasos_intermedios/matriz.npz", r"/data/pasos_intermedios/csv.csv")
npz_to_vcf(r"/data/pasos_intermedios/matriz.npz", r"/data/pasos_intermedios/vcf_transformado/vcf_transf.vcf")
split_vcf(r"/data/pasos_intermedios/vcf_transformado/vcf_transf.vcf", r"/data/pasos_intermedios/vcf_transformado/vcf_snv.vcf", r"/data/pasos_intermedios/vcf_transformado/vcf_sv.vcf")
