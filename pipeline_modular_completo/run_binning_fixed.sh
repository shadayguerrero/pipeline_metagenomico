#!/bin/bash

# Script para binning metagenómico con MetaBAT2, MaxBin2, CONCOCT y DAS Tool
# Versión completa con refinamiento de bins

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
MAGENTA='\033[0;35m'
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
    
    for tool in bowtie2-build bowtie2 samtools jgi_summarize_bam_contig_depths metabat2 run_MaxBin.pl cut_up_fasta.py concoct DAS_Tool; do
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
    mkdir -p "${SAMPLE_OUTPUT}"/{mapping,metabat2,maxbin2,concoct,dastool}
    cd "${SAMPLE_OUTPUT}"
    
    # ========================================================================
    # PASO 1: Filtrar contigs
    # ========================================================================
    print_step "[1/9] Filtrando contigs por longitud..."
    
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
    print_step "[2/9] Indexando contigs..."
    
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
    print_step "[3/9] Mapeando reads a contigs..."
    
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
    # PASO 4: Calcular profundidad para MetaBAT2
    # ========================================================================
    print_step "[4/9] Calculando profundidad de cobertura..."
    
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
    print_step "[5/9] Ejecutando binning con MetaBAT2..."
    
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
    # PASO 6: Binning con MaxBin2 (usando reads directamente)
    # ========================================================================
    print_step "[6/9] Ejecutando binning con MaxBin2..."
    
    # MaxBin2 usando reads directamente (más robusto que usar abundancia)
    run_MaxBin.pl \
        -contig "${FILTERED_CONTIGS}" \
        -reads "${R1_FILE}" \
        -reads2 "${R2_FILE}" \
        -out maxbin2/bin \
        -thread ${THREADS} \
        -min_contig_length ${MIN_BIN_SIZE} \
        > maxbin2/maxbin2.log 2>&1
    
    MAXBIN2_BINS=$(ls maxbin2/bin.*.fasta 2>/dev/null | wc -l)
    print_success "MaxBin2 completado: ${MAXBIN2_BINS} bins generados"
    
    # ========================================================================
    # PASO 8: Binning con CONCOCT
    # ========================================================================
    print_step "[8/9] Ejecutando binning con CONCOCT..."
    
    # Cortar contigs en fragmentos de 10kb
    cut_up_fasta.py \
        "${FILTERED_CONTIGS}" \
        -c 10000 \
        -o 0 \
        --merge_last \
        -b concoct/contigs_10K.bed \
        > concoct/contigs_10K.fa
    
    # Generar tabla de cobertura
    concoct_coverage_table.py \
        concoct/contigs_10K.bed \
        mapping/mapped_sorted.bam \
        > concoct/coverage_table.tsv
    
    # Ejecutar CONCOCT
    concoct \
        --composition_file concoct/contigs_10K.fa \
        --coverage_file concoct/coverage_table.tsv \
        -b concoct/ \
        -t ${THREADS} \
        > concoct/concoct.log 2>&1
    
    # Extraer bins
    merge_cutup_clustering.py \
        concoct/clustering_gt1000.csv \
        > concoct/clustering_merged.csv
    
    # Crear directorio de bins
    mkdir -p concoct/bins
    
    extract_fasta_bins.py \
        "${FILTERED_CONTIGS}" \
        concoct/clustering_merged.csv \
        --output_path concoct/bins/
    
    CONCOCT_BINS=$(ls concoct/bins/*.fa 2>/dev/null | wc -l)
    print_success "CONCOCT completado: ${CONCOCT_BINS} bins generados"
    
    # ========================================================================
    # PASO 9: Refinamiento con DAS Tool
    # ========================================================================
    print_step "[9/9] Refinando bins con DAS Tool..."
    
    # Crear archivos de scaffold-to-bin para cada binner
    # MetaBAT2
    if [ ${METABAT2_BINS} -gt 0 ]; then
        for bin in metabat2/bin.*.fa; do
            bin_name=$(basename "$bin" .fa)
            grep "^>" "$bin" | sed 's/>//' | awk -v bin="$bin_name" '{print $1"\t"bin}' >> dastool/metabat2_scaffolds2bin.tsv
        done
    fi
    
    # MaxBin2
    if [ ${MAXBIN2_BINS} -gt 0 ]; then
        for bin in maxbin2/bin.*.fasta; do
            bin_name=$(basename "$bin" .fasta)
            grep "^>" "$bin" | sed 's/>//' | awk -v bin="$bin_name" '{print $1"\t"bin}' >> dastool/maxbin2_scaffolds2bin.tsv
        done
    fi
    
    # CONCOCT
    if [ ${CONCOCT_BINS} -gt 0 ]; then
        for bin in concoct/bins/*.fa; do
            bin_name=$(basename "$bin" .fa)
            grep "^>" "$bin" | sed 's/>//' | awk -v bin="$bin_name" '{print $1"\t"bin}' >> dastool/concoct_scaffolds2bin.tsv
        done
    fi
    
    # Ejecutar DAS Tool
    DASTOOL_INPUT=""
    DASTOOL_LABELS=""
    
    [ -f dastool/metabat2_scaffolds2bin.tsv ] && DASTOOL_INPUT="${DASTOOL_INPUT},dastool/metabat2_scaffolds2bin.tsv" && DASTOOL_LABELS="${DASTOOL_LABELS},MetaBAT2"
    [ -f dastool/maxbin2_scaffolds2bin.tsv ] && DASTOOL_INPUT="${DASTOOL_INPUT},dastool/maxbin2_scaffolds2bin.tsv" && DASTOOL_LABELS="${DASTOOL_LABELS},MaxBin2"
    [ -f dastool/concoct_scaffolds2bin.tsv ] && DASTOOL_INPUT="${DASTOOL_INPUT},dastool/concoct_scaffolds2bin.tsv" && DASTOOL_LABELS="${DASTOOL_LABELS},CONCOCT"
    
    # Remover comas iniciales
    DASTOOL_INPUT=${DASTOOL_INPUT#,}
    DASTOOL_LABELS=${DASTOOL_LABELS#,}
    
    if [ -n "${DASTOOL_INPUT}" ]; then
        DAS_Tool \
            -i "${DASTOOL_INPUT}" \
            -l "${DASTOOL_LABELS}" \
            -c "${FILTERED_CONTIGS}" \
            -o dastool/DASTool \
            --write_bins \
            -t ${THREADS} \
            --search_engine diamond \
            > dastool/dastool.log 2>&1
        
        DASTOOL_BINS=$(ls dastool/DASTool_DASTool_bins/*.fa 2>/dev/null | wc -l)
        print_success "DAS Tool completado: ${DASTOOL_BINS} bins refinados"
    else
        print_warning "No hay bins para refinar con DAS Tool"
        DASTOOL_BINS=0
    fi
    
    # ========================================================================
    # RESUMEN
    # ========================================================================
    echo ""
    print_header "Resumen de ${SAMPLE}"
    echo -e "  ${CYAN}Contigs filtrados:${NC} ${FILTERED_COUNT}"
    echo -e "  ${CYAN}Reads mapeados:${NC} ${TOTAL_READS}"
    echo ""
    echo -e "  ${MAGENTA}Bins por binner:${NC}"
    echo -e "    ${GREEN}MetaBAT2:${NC} ${METABAT2_BINS} bins"
    echo -e "    ${GREEN}MaxBin2:${NC} ${MAXBIN2_BINS} bins"
    echo -e "    ${GREEN}CONCOCT:${NC} ${CONCOCT_BINS} bins"
    echo -e "    ${YELLOW}DAS Tool (refinados):${NC} ${DASTOOL_BINS} bins"
    echo ""
    echo -e "  ${CYAN}Directorios de salida:${NC}"
    echo -e "    metabat2/  - Bins de MetaBAT2"
    echo -e "    maxbin2/   - Bins de MaxBin2"
    echo -e "    concoct/   - Bins de CONCOCT"
    echo -e "    dastool/   - Bins refinados (RECOMENDADO)"
    echo ""
    
    return 0
}

# ============================================================================
# MAIN
# ============================================================================

print_header "BINNING METAGENÓMICO COMPLETO"
echo ""
echo "Configuración:"
echo "  Directorio de ensamblajes: ${ASSEMBLY_DIR}"
echo "  Directorio de reads: ${HOST_REMOVED_DIR}"
echo "  Directorio de salida: ${OUTPUT_DIR}"
echo "  Threads: ${THREADS}"
echo "  Longitud mínima de contigs: ${MIN_CONTIG_LEN} bp"
echo "  Tamaño mínimo de bins: ${MIN_BIN_SIZE} bp"
echo ""
echo "Binnners a ejecutar:"
echo "  1. MetaBAT2 - Binning basado en cobertura y composición"
echo "  2. MaxBin2  - Binning con marcadores de copia única"
echo "  3. CONCOCT  - Binning con fragmentación de contigs"
echo "  4. DAS Tool - Refinamiento y selección de mejores bins"
echo ""

# Verificar herramientas
print_step "Verificando herramientas necesarias..."
if ! check_tools; then
    echo ""
    print_error "Faltan herramientas necesarias"
    echo ""
    echo "Solución:"
    echo "  1. Activa el ambiente: micromamba activate binning"
    echo "  2. Verifica la instalación: micromamba list | grep -E 'metabat|maxbin|concoct|das_tool'"
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
echo -e "${YELLOW}Recomendación:${NC} Usa los bins de ${MAGENTA}dastool/DASTools_DASTool_bins/${NC}"
echo "para los análisis posteriores (GTDB-Tk, Prokka, etc.)"
echo ""

exit 0
