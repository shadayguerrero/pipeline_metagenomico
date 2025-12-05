#!/bin/bash

# ============================================================================
# Clasificación taxonómica con Kraken2 en paralelo
# Corrige problemas de delimitadores y variables de entorno
# ============================================================================

INPUT_DIR="/home_local/camda/shaday/4_cienegas/output/host_removed"
OUTPUT_DIR="/home_local/camda/shaday/4_cienegas/output/kraken2_results"
KRAKEN2_DB="/home_local/camda/shaday/database/k2_standard"

THREADS_PER_SAMPLE=8
MAX_PARALLEL_JOBS=1
TOTAL_THREADS=$((THREADS_PER_SAMPLE * MAX_PARALLEL_JOBS))

CONFIDENCE=0.0
MIN_BASE_QUALITY=0
MIN_HIT_GROUPS=2

SUFFIX_R1="_host_removed_R1.fastq.gz"
SUFFIX_R2="_host_removed_R2.fastq.gz"

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; CYAN='\033[0;36m'; NC='\033[0m'

process_kraken2() {
    local R1_FILE=$1
    local SAMPLE=$2
    local R2_FILE=$3
    local THREADS=$4
    local SAMPLE_OUTPUT="${OUTPUT_DIR}/${SAMPLE}"

    mkdir -p "${SAMPLE_OUTPUT}"

    echo -e "${BLUE}[$(date +%H:%M:%S)] Procesando muestra: ${YELLOW}${SAMPLE}${NC}"

    local KRAKEN_OUTPUT="${SAMPLE_OUTPUT}/${SAMPLE}.kraken"
    local REPORT_OUTPUT="${SAMPLE_OUTPUT}/${SAMPLE}.report"

    # Verificar que bowtie2 y kraken2 estén accesibles
    which kraken2 > "${SAMPLE_OUTPUT}/debug_path.txt" 2>&1

    # Ejecutar Kraken2
    kraken2 --db "${KRAKEN2_DB}" \
        --threads "${THREADS}" \
        --paired \
        --confidence "${CONFIDENCE}" \
        --minimum-base-quality "${MIN_BASE_QUALITY}" \
        --minimum-hit-groups "${MIN_HIT_GROUPS}" \
        --output "${KRAKEN_OUTPUT}" \
        --report "${REPORT_OUTPUT}" \
        --gzip-compressed \
        "${R1_FILE}" "${R2_FILE}" \
        2> "${SAMPLE_OUTPUT}/${SAMPLE}_kraken2.log"

    if [ $? -ne 0 ]; then
        echo -e "${RED}  ✗ Error en Kraken2 para ${SAMPLE}${NC}"
        return 1
    fi

    gzip -9 "${KRAKEN_OUTPUT}"

    echo -e "${GREEN}  ✓ Kraken2 completado para ${SAMPLE}${NC}"
    return 0
}

# Exportar variables necesarias
export -f process_kraken2
export OUTPUT_DIR KRAKEN2_DB THREADS_PER_SAMPLE CONFIDENCE MIN_BASE_QUALITY MIN_HIT_GROUPS RED GREEN YELLOW BLUE CYAN NC PATH

# Crear directorio de salida
mkdir -p "${OUTPUT_DIR}"

# Crear lista de muestras robusta
SAMPLE_LIST=$(mktemp)
for R1_FILE in "${INPUT_DIR}"/*/*"${SUFFIX_R1}"; do
    [ -f "$R1_FILE" ] || continue
    SAMPLE_DIR=$(dirname "$R1_FILE")
    SAMPLE=$(basename "$SAMPLE_DIR")
    R2_FILE="${R1_FILE/${SUFFIX_R1}/${SUFFIX_R2}}"
    if [ -f "$R2_FILE" ]; then
        echo "${R1_FILE}|${SAMPLE}|${R2_FILE}" >> "${SAMPLE_LIST}"
    else
        echo -e "${YELLOW}Advertencia: Falta ${R2_FILE}${NC}"
    fi
done

TOTAL_SAMPLES=$(wc -l < "${SAMPLE_LIST}")
[ "$TOTAL_SAMPLES" -eq 0 ] && { echo -e "${RED}No se encontraron muestras.${NC}"; exit 1; }

echo -e "${GREEN}Se procesarán ${TOTAL_SAMPLES} muestras.${NC}"
echo ""

# Procesamiento paralelo (corregido)
cat "${SAMPLE_LIST}" | parallel --colsep '\|' \
    --env PATH \
    --jobs "${MAX_PARALLEL_JOBS}" \
    --progress \
    --joblog "${OUTPUT_DIR}/parallel_joblog.txt" \
    process_kraken2 {1} {2} {3} "${THREADS_PER_SAMPLE}"

echo -e "${GREEN}Procesamiento completado.${NC}"
echo -e "Resultados: ${YELLOW}${OUTPUT_DIR}${NC}"

rm -f "${SAMPLE_LIST}"
