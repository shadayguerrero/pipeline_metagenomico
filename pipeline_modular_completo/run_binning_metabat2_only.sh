#!/bin/bash

# Script para binning metagenómico con MetaBAT2 únicamente
# Versión simplificada para ambientes sin MaxBin2/CONCOCT/DAS Tool

# ============================================================================
# CONFIGURACIÓN
# ============================================================================

# Detectar directorio base (donde está el script)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"

# Directorios (usa variables de entorno o valores por defecto)
ASSEMBLY_DIR="${ASSEMBLY_DIR:-${PROJECT_DIR}/output/megahit_assemblies}"
HOST_REMOVED_DIR="${HOST_REMOVED_DIR:-${PROJECT_DIR}/output/host_removed}"
OUTPUT_DIR="${OUTPUT_DIR:-${PROJECT_DIR}/output/binning}"

# Parámetros de procesamiento
THREADS="${THREADS:-8}"
MIN_CONTIG_LEN=500
MIN_BIN_SIZE=2500  # 2.5 kb para metagenomas

# Colores
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

# ============================================================================
# FUNCIONES
# ============================================================================

print_header() {
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}"
}

print_step() {
    echo -e "${CYAN}  $1${NC}"
}

print_success() {
    echo -e "${GREEN}    ✓ $1${NC}"
}

print_error() {
    echo -e "${RED}    ✗ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}    ⚠ $1${NC}"
}

# Verificar herramientas necesarias
check_tools() {
    local missing_tools=()
    
    for tool in bowtie2-build bowtie2 samtools jgi_summarize_bam_contig_depths metabat2; do
        if ! command -v $tool &> /dev/null; then
            missing_tools+=($tool)
        fi
    done
    
    if [ ${#missing_tools[@]} -gt 0 ]; then
        print_error "Herramientas faltantes: ${missing_tools[*]}"
        print_warning "Asegúrate de activar el ambiente 'binning'"
        print_warning "Comando: micromamba activate binning"
        return 1
    fi
    
    return 0
}

# Función para procesar una muestra
process_sample() {
    local SAMPLE=$1
    
    print_header "Procesando: ${SAMPLE}"
    
    local ASSEMBLY_PATH="${ASSEMBLY_DIR}/${SAMPLE}/final.contigs.fa"
    local R1_FILE="${HOST_REMOVED_DIR}/${SAMPLE}/${SAMPLE}_host_removed_R1.fastq.gz"
    local R2_FILE="${HOST_REMOVED_DIR}/${SAMPLE}/${SAMPLE}_host_removed_R2.fastq.gz"
    local SAMPLE_OUTPUT="${OUTPUT_DIR}/${SAMPLE}"
    
    # Verificar archivos
    if [ ! -f "${ASSEMBLY_PATH}" ]; then
        print_error "No se encontró el ensamblaje: ${ASSEMBLY_PATH}"
        return 1
    fi
    
    if [ ! -f "${R1_FILE}" ] || [ ! -f "${R2_FILE}" ]; then
        print_error "No se encontraron los reads"
        return 1
    fi
    
    # Crear directorios
    mkdir -p "${SAMPLE_OUTPUT}"/{mapping,metabat2}
    cd "${SAMPLE_OUTPUT}"
    
    # ========================================================================
    # PASO 1: Filtrar contigs
    # ========================================================================
    print_step "[1/5] Filtrando contigs por longitud..."
    
    FILTERED_CONTIGS="mapping/filtered_contigs.fa"
    
    awk -v min_len=${MIN_CONTIG_LEN} '
    BEGIN {write_seq = 0}
    /^>/ {
        if (match($0, /len=([0-9]+)/)) {
            contig_len = substr($0, RSTART+4, RLENGTH-4)
            write_seq = (contig_len >= min_len)
            if (write_seq) print
        }
        next
    }
    write_seq {print}
    ' "${ASSEMBLY_PATH}" > "${FILTERED_CONTIGS}"
    
    FILTERED_COUNT=$(grep -c "^>" "${FILTERED_CONTIGS}" 2>/dev/null || echo 0)
    print_success "${FILTERED_COUNT} contigs filtrados"
    
    if [ ${FILTERED_COUNT} -eq 0 ]; then
        print_error "No hay contigs >= ${MIN_CONTIG_LEN} bp"
        return 1
    fi
    
    # ========================================================================
    # PASO 2: Indexar contigs
    # ========================================================================
    print_step "[2/5] Indexando contigs..."
    
    bowtie2-build --threads ${THREADS} \
        "${FILTERED_CONTIGS}" \
        mapping/contigs_index \
        > mapping/bowtie2-build.log 2>&1
    
    if [ $? -ne 0 ]; then
        print_error "Error al indexar contigs"
        return 1
    fi
    
    print_success "Indexación completada"
    
    # ========================================================================
    # PASO 3: Mapear reads
    # ========================================================================
    print_step "[3/5] Mapeando reads a contigs..."
    
    bowtie2 \
        -x mapping/contigs_index \
        -1 "${R1_FILE}" \
        -2 "${R2_FILE}" \
        -p ${THREADS} \
        --no-unal \
        2> mapping/bowtie2.log | \
    samtools sort -@ ${THREADS} -o mapping/mapped_sorted.bam -
    
    if [ $? -ne 0 ]; then
        print_error "Error al mapear reads"
        return 1
    fi
    
    samtools index -@ ${THREADS} mapping/mapped_sorted.bam
    
    TOTAL_READS=$(samtools view -c mapping/mapped_sorted.bam)
    print_success "Mapeo completado: ${TOTAL_READS} reads mapeados"
    
    # ========================================================================
    # PASO 4: Calcular profundidad
    # ========================================================================
    print_step "[4/5] Calculando profundidad de cobertura..."
    
    jgi_summarize_bam_contig_depths \
        --outputDepth metabat2/depth.txt \
        mapping/mapped_sorted.bam \
        > metabat2/depth.log 2>&1
    
    if [ $? -ne 0 ]; then
        print_error "Error al calcular profundidad"
        return 1
    fi
    
    print_success "Profundidad calculada"
    
    # ========================================================================
    # PASO 5: Binning con MetaBAT2
    # ========================================================================
    print_step "[5/5] Ejecutando binning con MetaBAT2..."
    
    metabat2 \
        -i "${FILTERED_CONTIGS}" \
        -a metabat2/depth.txt \
        -o metabat2/bin \
        -m ${MIN_BIN_SIZE} \
        -t ${THREADS} \
        --seed 42 \
        > metabat2/metabat2.log 2>&1
    
    METABAT2_BINS=$(ls metabat2/bin.*.fa 2>/dev/null | wc -l)
    print_success "MetaBAT2 completado: ${METABAT2_BINS} bins generados"
    
    # ========================================================================
    # RESUMEN
    # ========================================================================
    echo ""
    print_header "Resumen de ${SAMPLE}"
    echo -e "  ${CYAN}Contigs filtrados:${NC} ${FILTERED_COUNT}"
    echo -e "  ${CYAN}Reads mapeados:${NC} ${TOTAL_READS}"
    echo -e "  ${GREEN}Bins generados:${NC} ${METABAT2_BINS}"
    echo ""
    echo -e "  ${CYAN}Directorio de salida:${NC}"
    echo -e "    metabat2/  - Bins de MetaBAT2"
    echo ""
    
    return 0
}

# ============================================================================
# MAIN
# ============================================================================

print_header "BINNING METAGENÓMICO - MetaBAT2"
echo ""
echo "Configuración:"
echo "  Directorio de ensamblajes: ${ASSEMBLY_DIR}"
echo "  Directorio de reads: ${HOST_REMOVED_DIR}"
echo "  Directorio de salida: ${OUTPUT_DIR}"
echo "  Threads: ${THREADS}"
echo "  Longitud mínima de contigs: ${MIN_CONTIG_LEN} bp"
echo "  Tamaño mínimo de bins: ${MIN_BIN_SIZE} bp"
echo ""
echo "Binner a ejecutar:"
echo "  - MetaBAT2 (binning basado en cobertura y composición)"
echo ""
echo -e "${YELLOW}Nota:${NC} Esta es la versión simplificada que solo usa MetaBAT2."
echo "Para usar MaxBin2, CONCOCT y DAS Tool, ejecuta: bash actualizar_ambiente_binning.sh"
echo ""

# Verificar herramientas
print_step "Verificando herramientas necesarias..."
if ! check_tools; then
    echo ""
    print_error "Faltan herramientas necesarias"
    echo ""
    echo "Solución:"
    echo "  1. Activa el ambiente: micromamba activate binning"
    echo "  2. Verifica la instalación: micromamba list | grep -E 'bowtie2|samtools|metabat'"
    echo ""
    exit 1
fi
print_success "Todas las herramientas disponibles"
echo ""

# Verificar directorios
if [ ! -d "${ASSEMBLY_DIR}" ]; then
    print_error "No se encontró el directorio de ensamblajes: ${ASSEMBLY_DIR}"
    exit 1
fi

if [ ! -d "${HOST_REMOVED_DIR}" ]; then
    print_error "No se encontró el directorio de reads: ${HOST_REMOVED_DIR}"
    exit 1
fi

# Crear directorio de salida
mkdir -p "${OUTPUT_DIR}"

# Encontrar muestras
SAMPLES=$(ls -d ${ASSEMBLY_DIR}/*/ 2>/dev/null | xargs -n 1 basename)

if [ -z "${SAMPLES}" ]; then
    print_error "No se encontraron muestras en ${ASSEMBLY_DIR}"
    exit 1
fi

echo "Muestras encontradas:"
for sample in ${SAMPLES}; do
    echo "  - ${sample}"
done
echo ""

# Procesar cada muestra
SUCCESS_COUNT=0
FAIL_COUNT=0

for sample in ${SAMPLES}; do
    if process_sample "${sample}"; then
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    else
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
    echo ""
done

# Resumen final
print_header "BINNING COMPLETADO"
echo ""
echo "Resultados:"
echo -e "  ${GREEN}✓ Exitosas: ${SUCCESS_COUNT}${NC}"
echo -e "  ${RED}✗ Fallidas: ${FAIL_COUNT}${NC}"
echo ""
echo "Los bins están en: ${OUTPUT_DIR}"
echo ""
echo -e "${YELLOW}Siguiente paso:${NC}"
echo "  Ejecuta el módulo 6 (GTDB-Tk) para clasificar los bins:"
echo "  bash run_gtdbtk.sh"
echo ""

exit 0
