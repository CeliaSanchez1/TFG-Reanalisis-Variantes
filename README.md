# Trabajo de Fin de Grado: Implementación de un sistema automatizado de reanotación periódica de exomas completos para mejorar la tasa diagnóstica en un entorno asistencial.

# Descripción
Este repositorio contiene el pipeline desarrollado para el análisis de variantes genómicas a partir de datos de exoma. Se recogen scripts para fusionar y anotar archivos VCF, así como para comprobar y descargar las últimas actualizaciones de las bases de datos y comparar anotaciones con distintas versiones de estas. 
El repositorio incluye tanto el código del pipeline (en distintas carpetas para los scripts en lenguaje Python y Bash) como los resultados generados en diferentes etapas del análisis (en otra carpeta aislada de resultados). Cabe recalcar que estos resultados son los obtenidos al trabajar sobre archivos VCF y versiones de las bases de datos concretas, serán diferentes según se definan estos factores en la ejecución.

# Flujo de trabajo
El pipeline sigue las siguientes etapas:

1. Comprobación de versiones de las bases de datos
2. Integración del nuevo lote en el msVCF
3. División del msVCF en variantes cortas y largas
4. División del msVCF por cromosomas
5. Reanotación de los cromosomas
6. Reanotación del cromosoma M
7. Comparación de versiones de la anotación de los cromosomas
8. Comparación de versiones de la anotación del cromosoma M
9. Detección de nuevas variantes patogénicas
10. Filtrado por HPO y extracción de TSV con campos más relevantes
11. Unión de los TSV de todos los cromosomas
12. Integración con información diagnóstica

# Ejecución
Todo el flujo de trabajo se gestiona a través del script workflow.sh

# Requisitos
Todas las librerías empleadas se recogen en un entorno gestionado con Conda. Creación y activación:

```bash
conda env create -f environment.yml
conda activate tfg-env
```

# Autor
Celia Sánchez Tejeda – Desarrollado en Universidad Pública de Navarra y Navarra de Servicios y Tecnologías S. A.
