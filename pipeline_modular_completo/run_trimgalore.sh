#!/bin/bash

# Script para ejecutar Trim Galore en archivos FASTQ pareados
# Autor: Script generado para análisis de 4 Ciénegas
# Fecha: $(date +%Y-%m-%d)

# Configuración de directorios
INPUT_DIR="/home_local/camda/shaday/4_cienegas/MergedFastq"
OUTPUT_DIR="/home_local/camda/shaday/4_cienegas/output"

# Crear directorio de salida si no existe
mkdir -p "${OUTPUT_DIR}"

# Parámetros de Trim Galore
THREADS=20  # Ajustar según recursos disponibles
QUALITY=20  # Calidad mínima de Phred score
LENGTH=20   # Longitud mínima de reads después del trimming

# Colores para output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  Trim Galore - Procesamiento Pareado  ${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo -e "Directorio de entrada: ${YELLOW}${INPUT_DIR}${NC}"
echo -e "Directorio de salida: ${YELLOW}${OUTPUT_DIR}${NC}"
echo -e "Threads: ${YELLOW}${THREADS}${NC}"
echo -e "Calidad mínima: ${YELLOW}${QUALITY}${NC}"
echo -e "Longitud mínima: ${YELLOW}${LENGTH}${NC}"
echo ""

# Verificar que existe el directorio de entrada
if [ ! -d "${INPUT_DIR}" ]; then
    echo -e "${RED}Error: El directorio de entrada no existe: ${INPUT_DIR}${NC}"
    exit 1
fi

# Verificar que Trim Galore está instalado
if ! command -v trim_galore &> /dev/null; then
    echo -e "${RED}Error: Trim Galore no está instalado${NC}"
    echo "Instalar con: conda install -c bioconda trim-galore"
    exit 1
fi

# Contador de muestras
TOTAL=0
SUCCESS=0
FAILED=0

# Procesar archivos pareados
# Identificar muestras únicas (sin sufijos _1 y _2)
cd "${INPUT_DIR}"

for R1_FILE in *_1.fastq.gz; do
    # Extraer nombre base de la muestra
    SAMPLE=$(basename "${R1_FILE}" _1.fastq.gz)
    R2_FILE="${SAMPLE}_2.fastq.gz"
    
    # Verificar que existe el archivo R2 correspondiente
    if [ ! -f "${R2_FILE}" ]; then
        echo -e "${RED}Advertencia: No se encontró el archivo pareado ${R2_FILE} para ${R1_FILE}${NC}"
        ((FAILED++))
        continue
    fi
    
    ((TOTAL++))
    
    echo -e "${GREEN}[${TOTAL}] Procesando muestra: ${YELLOW}${SAMPLE}${NC}"
    echo -e "    R1: ${R1_FILE}"
    echo -e "    R2: ${R2_FILE}"
    
    # Ejecutar Trim Galore en modo pareado
    trim_galore \
        --paired \
        --quality ${QUALITY} \
        --length ${LENGTH} \
        --cores ${THREADS} \
        --fastqc \
        --output_dir "${OUTPUT_DIR}" \
        "${INPUT_DIR}/${R1_FILE}" \
        "${INPUT_DIR}/${R2_FILE}"
    
    # Verificar si el comando fue exitoso
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}    ✓ Completado exitosamente${NC}"
        ((SUCCESS++))
    else
        echo -e "${RED}    ✗ Error en el procesamiento${NC}"
        ((FAILED++))
    fi
    
    echo ""
done

# Resumen final
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  Resumen del Procesamiento${NC}"
echo -e "${GREEN}========================================${NC}"
echo -e "Total de muestras procesadas: ${YELLOW}${TOTAL}${NC}"
echo -e "Exitosas: ${GREEN}${SUCCESS}${NC}"
echo -e "Fallidas: ${RED}${FAILED}${NC}"
echo ""
echo -e "${GREEN}Los archivos procesados se encuentran en:${NC}"
echo -e "${YELLOW}${OUTPUT_DIR}${NC}"
echo ""
echo -e "${GREEN}Proceso completado$(date +%Y-%m-%d\ %H:%M:%S)${NC}"
