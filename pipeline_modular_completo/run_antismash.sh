#!/bin/bash

# ==============================================================
# AntiSMASH: Análisis de clusters de genes de metabolitos secundarios
# Requiere: Prokka completado (usa archivos .gbk)
# ==============================================================

# Detectar directorio base (donde está el script)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"

# Directorios (usa variables de entorno o valores por defecto)
PROKKA_DIR="${PROKKA_DIR:-${PROJECT_DIR}/output/prokka}"
OUTPUT_DIR="${OUTPUT_DIR:-${PROJECT_DIR}/output/antismash}"
THREADS="${THREADS:-8}"

# Selección de bins (debe coincidir con Prokka)
# Opciones: dastool, metabat2, maxbin2, concoct, auto
BINS_SOURCE="${BINS_SOURCE:-auto}"

# Parámetros de AntiSMASH
GENEFINDING_TOOL="none"  # Usa anotaciones de Prokka
MIN_CONTIG_LEN=1000      # Longitud mínima de contigs para analizar

# Colores
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

# ==============================================================
# FUNCIONES
# ==============================================================

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

# ==============================================================
# PREPARACIÓN
# ==============================================================

print_header "AntiSMASH - Metabolitos Secundarios"
echo ""
echo "Configuración:"
echo "  Directorio Prokka: ${PROKKA_DIR}"
echo "  Directorio de salida: ${OUTPUT_DIR}"
echo "  Fuente de bins: ${BINS_SOURCE}"
echo "  Threads: ${THREADS}"
echo "  Longitud mínima de contigs: ${MIN_CONTIG_LEN} bp"
echo ""

# Verificar que AntiSMASH esté disponible
if ! command -v antismash &> /dev/null; then
    print_error "AntiSMASH no está disponible"
    echo ""
    echo "Solución:"
    echo "  1. Activa el ambiente: micromamba activate antismash"
    echo "  2. Verifica la instalación: antismash --version"
    echo "  3. Descarga las bases de datos: download-antismash-databases"
    echo ""
    exit 1
fi

ANTISMASH_VERSION=$(antismash --version 2>&1 | head -n1)
print_success "AntiSMASH detectado: ${ANTISMASH_VERSION}"
echo ""

# Crear directorio de salida
mkdir -p "${OUTPUT_DIR}"

# Verificar que Prokka se haya ejecutado
if [ ! -d "${PROKKA_DIR}" ]; then
    print_error "No se encontró el directorio de Prokka: ${PROKKA_DIR}"
    echo ""
    echo "Ejecuta primero: bash run_prokka.sh"
    echo ""
    exit 1
fi

# Encontrar muestras
SAMPLES=$(ls -d ${PROKKA_DIR}/*/ 2>/dev/null | xargs -n 1 basename)

if [ -z "${SAMPLES}" ]; then
    print_error "No se encontraron muestras en ${PROKKA_DIR}"
    exit 1
fi

echo "Muestras encontradas:"
for sample in ${SAMPLES}; do
    echo "  - ${sample}"
done
echo ""

# ==============================================================
# PROCESAMIENTO
# ==============================================================

SUCCESS_COUNT=0
FAIL_COUNT=0
TOTAL_BINS=0
TOTAL_CLUSTERS=0

for SAMPLE in ${SAMPLES}; do
    print_header "Procesando: ${SAMPLE}"
    
    SAMPLE_PROKKA="${PROKKA_DIR}/${SAMPLE}"
    SAMPLE_OUT="${OUTPUT_DIR}/${SAMPLE}"
    mkdir -p "${SAMPLE_OUT}"
    
    # Encontrar bins anotados
    BINS=$(ls -d ${SAMPLE_PROKKA}/*/ 2>/dev/null | xargs -n 1 basename)
    
    if [ -z "${BINS}" ]; then
        print_warning "No se encontraron bins anotados para ${SAMPLE}"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        echo ""
        continue
    fi
    
    BIN_COUNT=$(echo "${BINS}" | wc -l)
    print_step "Bins anotados: ${BIN_COUNT}"
    echo ""
    
    SAMPLE_SUCCESS=0
    SAMPLE_CLUSTERS=0
    
    for BIN in ${BINS}; do
        BIN_PROKKA="${SAMPLE_PROKKA}/${BIN}"
        BIN_OUT="${SAMPLE_OUT}/${BIN}"
        
        # Buscar archivo GenBank de Prokka
        GBK_FILE="${BIN_PROKKA}/${BIN}.gbk"
        
        if [ ! -f "${GBK_FILE}" ]; then
            print_warning "${BIN}: No se encontró archivo GenBank (.gbk)"
            continue
        fi
        
        print_step "Analizando: ${BIN}"
        
        # Ejecutar AntiSMASH
        antismash \
            --genefinding-tool "${GENEFINDING_TOOL}" \
            --output-dir "${BIN_OUT}" \
            --cpus ${THREADS} \
            --minlength ${MIN_CONTIG_LEN} \
            --skip-zip-file \
            "${GBK_FILE}" \
            > "${BIN_OUT}/antismash.log" 2>&1
        
        if [ $? -eq 0 ]; then
            # Contar clusters
            if [ -f "${BIN_OUT}/index.html" ]; then
                # Buscar clusters en el archivo JSON
                JSON_FILE=$(find "${BIN_OUT}" -name "*.json" | head -n1)
                if [ -f "${JSON_FILE}" ]; then
                    CLUSTERS=$(grep -o '"type":' "${JSON_FILE}" | wc -l)
                    if [ ${CLUSTERS} -gt 0 ]; then
                        print_success "${CLUSTERS} clusters de BGC encontrados"
                        SAMPLE_CLUSTERS=$((SAMPLE_CLUSTERS + CLUSTERS))
                    else
                        print_success "Sin clusters de BGC"
                    fi
                else
                    print_success "Análisis completado"
                fi
                SAMPLE_SUCCESS=$((SAMPLE_SUCCESS + 1))
                TOTAL_BINS=$((TOTAL_BINS + 1))
            else
                print_warning "AntiSMASH completado pero sin index.html"
            fi
        else
            print_error "Error en AntiSMASH (ver log)"
        fi
    done
    
    if [ ${SAMPLE_SUCCESS} -gt 0 ]; then
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
        TOTAL_CLUSTERS=$((TOTAL_CLUSTERS + SAMPLE_CLUSTERS))
        echo ""
        print_success "${SAMPLE}: ${SAMPLE_SUCCESS} bins analizados, ${SAMPLE_CLUSTERS} clusters BGC"
    else
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
    
    echo ""
done

# ==============================================================
# GENERAR RESUMEN
# ==============================================================

print_header "Generando Resumen"
echo ""

SUMMARY="${OUTPUT_DIR}/AntiSMASH_Summary.txt"

{
    echo "AntiSMASH - Resumen de Clusters de Genes Biosintéticos"
    echo "======================================================="
    echo ""
    echo "Muestras analizadas: ${SUCCESS_COUNT}"
    echo "Bins analizados: ${TOTAL_BINS}"
    echo "Clusters BGC totales: ${TOTAL_CLUSTERS}"
    echo ""
    echo "Resultados por muestra:"
    echo ""
    
    for SAMPLE in ${SAMPLES}; do
        SAMPLE_OUT="${OUTPUT_DIR}/${SAMPLE}"
        if [ -d "${SAMPLE_OUT}" ]; then
            BINS=$(ls -d ${SAMPLE_OUT}/*/ 2>/dev/null | xargs -n 1 basename)
            BIN_COUNT=$(echo "${BINS}" | wc -l)
            
            SAMPLE_CLUSTERS=0
            for BIN in ${BINS}; do
                JSON_FILE=$(find "${SAMPLE_OUT}/${BIN}" -name "*.json" 2>/dev/null | head -n1)
                if [ -f "${JSON_FILE}" ]; then
                    CLUSTERS=$(grep -o '"type":' "${JSON_FILE}" | wc -l)
                    SAMPLE_CLUSTERS=$((SAMPLE_CLUSTERS + CLUSTERS))
                fi
            done
            
            echo "  ${SAMPLE}: ${BIN_COUNT} bins, ${SAMPLE_CLUSTERS} clusters"
        fi
    done
    
    echo ""
    echo "Visualización:"
    echo "  Abre los archivos index.html en cada directorio de bin"
    echo "  Ejemplo: ${OUTPUT_DIR}/SAMPLE/BIN/index.html"
    echo ""
} > "${SUMMARY}"

cat "${SUMMARY}"

# ==============================================================
# RESUMEN FINAL
# ==============================================================

print_header "ANTISMASH COMPLETADO"
echo ""
echo "Resultados:"
echo -e "  ${GREEN}✓ Muestras exitosas: ${SUCCESS_COUNT}${NC}"
echo -e "  ${RED}✗ Muestras fallidas: ${FAIL_COUNT}${NC}"
echo -e "  ${CYAN}Total de bins analizados: ${TOTAL_BINS}${NC}"
echo -e "  ${YELLOW}Clusters BGC: ${TOTAL_CLUSTERS}${NC}"
echo ""
echo "Resultados individuales en: ${OUTPUT_DIR}"
echo "Resumen: ${SUMMARY}"
echo ""
echo -e "${CYAN}Para visualizar los resultados:${NC}"
echo "  Abre los archivos HTML en un navegador:"
echo "  firefox ${OUTPUT_DIR}/SAMPLE/BIN/index.html"
echo ""

exit 0
