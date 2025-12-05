#!/bin/bash

# ============================================================================
# Script para remover reads del host usando Bowtie2 y GNU Parallel
# Procesa archivos FASTQ pareados del output de Trim Galore
# Autor: Shaday Guerrero (versión corregida)
# Fecha: $(date +%Y-%m-%d)
# ============================================================================

# ============================================================================ #
# CONFIGURACIÓN
# ============================================================================ #

INPUT_DIR="/home_local/camda/shaday/4_cienegas/output/trim"        # Directorio de entrada
OUTPUT_DIR="/home_local/camda/shaday/4_cienegas/output/host_removed"  # Directorio de salida
BOWTIE2_INDEX="/home_local/camda/shaday/database/bowtie2/index/All/all"  # Índice de Bowtie2

THREADS_PER_SAMPLE=12   # Threads por muestra individual
MAX_PARALLEL_JOBS=2     # Número de muestras a procesar en paralelo
TOTAL_THREADS=$((THREADS_PER_SAMPLE * MAX_PARALLEL_JOBS))

SUFFIX_R1="_1_val_1.fq.gz"
SUFFIX_R2="_2_val_2.fq.gz"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# ============================================================================ #
# FUNCIÓN PRINCIPAL PARA PROCESAR UNA MUESTRA
# ============================================================================ #

process_sample() {
    local R1_FILE=$1
    local SAMPLE=$2
    local R2_FILE=$3
    local THREADS=$4
    local SAMPLE_OUTPUT="${OUTPUT_DIR}/${SAMPLE}"

    mkdir -p "${SAMPLE_OUTPUT}"

    echo -e "${BLUE}[$(date +%H:%M:%S)] Procesando muestra: ${YELLOW}${SAMPLE}${NC}"

    # Verificar disponibilidad de Bowtie2 y Samtools
    which bowtie2 > "${SAMPLE_OUTPUT}/debug_path.txt" 2>&1
    which samtools >> "${SAMPLE_OUTPUT}/debug_path.txt" 2>&1

    # PASO 1: Mapeo
    echo -e "${BLUE}  [1/5] Mapeando reads contra el host...${NC}"
    bowtie2 -p "${THREADS}" -x "${BOWTIE2_INDEX}" \
        -1 "${R1_FILE}" -2 "${R2_FILE}" \
        -S "${SAMPLE_OUTPUT}/${SAMPLE}_results.sam" \
        2> "${SAMPLE_OUTPUT}/bowtie2_verbose.txt"

    if [ $? -ne 0 ]; then
        echo -e "${RED}  ✗ Error en mapeo de ${SAMPLE}${NC}"
        return 1
    fi

    # PASO 2: SAM → BAM
    echo -e "${BLUE}  [2/5] Convirtiendo SAM a BAM...${NC}"
    samtools view -bS "${SAMPLE_OUTPUT}/${SAMPLE}_results.sam" \
        --threads "${THREADS}" > "${SAMPLE_OUTPUT}/${SAMPLE}_results.bam"
    if [ $? -ne 0 ]; then
        echo -e "${RED}  ✗ Error en conversión SAM→BAM de ${SAMPLE}${NC}"
        return 1
    fi

    # PASO 3: Filtrar reads no mapeados
    echo -e "${BLUE}  [3/5] Filtrando reads no mapeados (sin host)...${NC}"
    samtools view -b -f 13 -F 256 \
        "${SAMPLE_OUTPUT}/${SAMPLE}_results.bam" \
        --threads "${THREADS}" > "${SAMPLE_OUTPUT}/${SAMPLE}_bothEndsUnmapped.bam"
    if [ $? -ne 0 ]; then
        echo -e "${RED}  ✗ Error en filtrado de ${SAMPLE}${NC}"
        return 1
    fi

    # PASO 4: Ordenar BAM por nombre
    echo -e "${BLUE}  [4/5] Ordenando BAM por nombre...${NC}"
    samtools sort -n -@ "${THREADS}" \
        "${SAMPLE_OUTPUT}/${SAMPLE}_bothEndsUnmapped.bam" \
        -o "${SAMPLE_OUTPUT}/${SAMPLE}_sorted.bam"
    if [ $? -ne 0 ]; then
        echo -e "${RED}  ✗ Error en ordenamiento de ${SAMPLE}${NC}"
        return 1
    fi

    # PASO 5: Convertir a FASTQ
    echo -e "${BLUE}  [5/5] Generando FASTQ sin host...${NC}"
    samtools fastq -@ "${THREADS}" \
        "${SAMPLE_OUTPUT}/${SAMPLE}_sorted.bam" \
        -1 "${SAMPLE_OUTPUT}/${SAMPLE}_host_removed_R1.fastq.gz" \
        -2 "${SAMPLE_OUTPUT}/${SAMPLE}_host_removed_R2.fastq.gz"
    if [ $? -ne 0 ]; then
        echo -e "${RED}  ✗ Error en conversión a FASTQ de ${SAMPLE}${NC}"
        return 1
    fi

    # LIMPIEZA
    rm -f "${SAMPLE_OUTPUT}/${SAMPLE}_results.sam" \
          "${SAMPLE_OUTPUT}/${SAMPLE}_results.bam" \
          "${SAMPLE_OUTPUT}/${SAMPLE}_bothEndsUnmapped.bam"
    gzip -9 "${SAMPLE_OUTPUT}/${SAMPLE}_sorted.bam"

    # ESTADÍSTICAS
    echo -e "${BLUE}  [Estadísticas] Generando resumen...${NC}"
    R1_COUNT=$(zcat "${SAMPLE_OUTPUT}/${SAMPLE}_host_removed_R1.fastq.gz" | wc -l)
    R1_READS=$((R1_COUNT / 4))
    TOTAL_READS=$(grep "reads; of these:" "${SAMPLE_OUTPUT}/bowtie2_verbose.txt" | awk '{print $1}')
    OVERALL_ALIGNMENT=$(grep "overall alignment rate" "${SAMPLE_OUTPUT}/bowtie2_verbose.txt" | awk '{print $1}')

    cat > "${SAMPLE_OUTPUT}/summary_stats.txt" << EOF
===========================================
Resumen de Remoción de Host - ${SAMPLE}
===========================================
Fecha: $(date +%Y-%m-%d\ %H:%M:%S)
- Total de pares de reads procesados: ${TOTAL_READS}
- Tasa de alineamiento general: ${OVERALL_ALIGNMENT}
- Reads sin host (R1): ${R1_READS}
EOF

    echo -e "${GREEN}  ✓ ${SAMPLE} completado (${R1_READS} reads sin host)${NC}"
    return 0
}

# ============================================================================ #
# EXPORTAR VARIABLES Y FUNCIONES PARA GNU PARALLEL
# ============================================================================ #
export -f process_sample
export INPUT_DIR OUTPUT_DIR BOWTIE2_INDEX THREADS_PER_SAMPLE RED GREEN YELLOW BLUE NC PATH

# ============================================================================ #
# VERIFICACIONES PREVIAS
# ============================================================================ #

if [ ! -d "${INPUT_DIR}" ]; then
    echo -e "${RED}Error: Directorio de entrada no existe: ${INPUT_DIR}${NC}"
    exit 1
fi

if [ ! -f "${BOWTIE2_INDEX}.1.bt2l" ] && [ ! -f "${BOWTIE2_INDEX}.1.bt2" ]; then
    echo -e "${RED}Error: No se encontró el índice de Bowtie2: ${BOWTIE2_INDEX}${NC}"
    exit 1
fi

command -v bowtie2 >/dev/null 2>&1 || { echo -e "${RED}Error: bowtie2 no instalado${NC}"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo -e "${RED}Error: samtools no instalado${NC}"; exit 1; }

if ! command -v parallel &> /dev/null; then
    echo -e "${YELLOW}Advertencia: GNU Parallel no instalado. Se usará modo secuencial.${NC}"
    USE_PARALLEL=false
else
    USE_PARALLEL=true
fi

mkdir -p "${OUTPUT_DIR}"

# ============================================================================ #
# CREAR LISTA DE MUESTRAS
# ============================================================================ #

SAMPLE_LIST=$(mktemp)
cd "${INPUT_DIR}"

for R1_FILE in *${SUFFIX_R1}; do
    [ -f "${R1_FILE}" ] || continue
    SAMPLE=$(basename "${R1_FILE}" ${SUFFIX_R1})
    R2_FILE="${SAMPLE}${SUFFIX_R2}"
    [ -f "${R2_FILE}" ] || { echo -e "${RED}Advertencia: Falta ${R2_FILE}${NC}"; continue; }
    echo "${INPUT_DIR}/${R1_FILE}|${SAMPLE}|${INPUT_DIR}/${R2_FILE}" >> "${SAMPLE_LIST}"
done

TOTAL_SAMPLES=$(wc -l < "${SAMPLE_LIST}")
[ "${TOTAL_SAMPLES}" -eq 0 ] && { echo -e "${RED}No hay muestras.${NC}"; exit 1; }

echo -e "${GREEN}Se procesarán ${TOTAL_SAMPLES} muestras.${NC}"

# ============================================================================ #
# PROCESAMIENTO
# ============================================================================ #

START_TIME=$(date +%s)

if [ "${USE_PARALLEL}" = true ]; then
    echo -e "${GREEN}Iniciando procesamiento en paralelo...${NC}"
    cat "${SAMPLE_LIST}" | parallel --colsep '\|' \
        --env PATH \
        --jobs ${MAX_PARALLEL_JOBS} \
        --progress \
        --joblog "${OUTPUT_DIR}/parallel_joblog.txt" \
        process_sample {1} {2} {3} ${THREADS_PER_SAMPLE}
else
    COUNTER=0
    while IFS='|' read -r R1_FILE SAMPLE R2_FILE; do
        ((COUNTER++))
        echo -e "${GREEN}[${COUNTER}/${TOTAL_SAMPLES}]${NC}"
        process_sample "${R1_FILE}" "${SAMPLE}" "${R2_FILE}" ${THREADS_PER_SAMPLE}
    done < "${SAMPLE_LIST}"
fi

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
ELAPSED_MIN=$((ELAPSED / 60))
ELAPSED_SEC=$((ELAPSED % 60))

# ============================================================================ #
# RESUMEN FINAL
# ============================================================================ #

SUMMARY_FILE="${OUTPUT_DIR}/00_RESUMEN_GENERAL.txt"
echo "=========================================" > "${SUMMARY_FILE}"
echo "RESUMEN GENERAL - REMOCIÓN DE HOST" >> "${SUMMARY_FILE}"
echo "=========================================" >> "${SUMMARY_FILE}"
echo "Fecha: $(date +%Y-%m-%d\ %H:%M:%S)" >> "${SUMMARY_FILE}"
echo "Total de muestras: ${TOTAL_SAMPLES}" >> "${SUMMARY_FILE}"
echo "Tiempo total: ${ELAPSED_MIN}m ${ELAPSED_SEC}s" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

for SAMPLE_DIR in "${OUTPUT_DIR}"/*/; do
    [ -d "${SAMPLE_DIR}" ] || continue
    SAMPLE=$(basename "${SAMPLE_DIR}")
    if [ -f "${SAMPLE_DIR}/summary_stats.txt" ]; then
        READS=$(grep "Reads sin host" "${SAMPLE_DIR}/summary_stats.txt" | awk '{print $5}')
        RATE=$(grep "Tasa de alineamiento" "${SAMPLE_DIR}/summary_stats.txt" | cut -d':' -f2 | xargs)
        echo "${SAMPLE}: ${READS} reads sin host (Alineamiento:${RATE})" >> "${SUMMARY_FILE}"
    fi
done

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  Proceso completado${NC}"
echo -e "${GREEN}  Total: ${TOTAL_SAMPLES} muestras en ${ELAPSED_MIN}m ${ELAPSED_SEC}s${NC}"
echo -e "${GREEN}  Resumen: ${SUMMARY_FILE}${NC}"

rm -f "${SAMPLE_LIST}"
exit 0

