# Trabajo de Fin de Grado- Reanotación de variantes genómicas

# Descripción
Este repositorio contiene el pipeline desarrollado para el análisis de variantes genómicas a partir de datos de exoma. Se recogen scripts para fusionar y anotar archivos VCF, así como para comprobar y descargar las últimas actualizaciones de las bases de datos y comparar anotaciones con distintas versiones de estas. 
El repositorio incluye tanto el código del pipeline como los resultados generados en diferentes etapas del análisis.

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
9. Unión de los archivos de todos los cromosomas
10. Detección de nuevas variantes patogénicas
11. Filtrar por HPO e integrar información diagnóstica

# Ejecución
Todo el flujo de trabajo se gestiona a través del script workflow.sh

# Requisitos
pip install -r requirements.txt

# Autor
Celia Sánchez Tejeda – Desarrollado en Universidad Pública de Navarra y Navarra de Servicios y Tecnologías S. A.
