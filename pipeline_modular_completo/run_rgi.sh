#!/bin/bash

# ==============================================================
# RGI: Análisis de resistencia a antibióticos usando CARD
# Requiere: Prokka completado
# ==============================================================

# Detectar directorio base (donde está el script)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"

# Directorios (usa variables de entorno o valores por defecto)
PROKKA_DIR="${PROKKA_DIR:-${PROJECT_DIR}/output/prokka}"
OUTPUT_DIR="${OUTPUT_DIR:-${PROJECT_DIR}/output/rgi}"
THREADS="${THREADS:-8}"

# Selección de bins (debe coincidir con Prokka)
# Opciones: dastool, metabat2, maxbin2, concoct, auto
BINS_SOURCE="${BINS_SOURCE:-auto}"

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

print_header "RGI - Resistencia a Antibióticos (CARD)"
echo ""
echo "Configuración:"
echo "  Directorio Prokka: ${PROKKA_DIR}"
echo "  Directorio de salida: ${OUTPUT_DIR}"
echo "  Fuente de bins: ${BINS_SOURCE}"
echo "  Threads: ${THREADS}"
echo ""

# Verificar que RGI esté disponible
if ! command -v rgi &> /dev/null; then
    print_error "RGI no está disponible"
    echo ""
    echo "Solución:"
    echo "  1. Activa el ambiente: micromamba activate rgi"
    echo "  2. Verifica la instalación: rgi main --version"
    echo "  3. Carga la base de datos CARD: rgi load --card_json /ruta/a/card.json"
    echo ""
    exit 1
fi

RGI_VERSION=$(rgi main --version 2>&1 | head -n1)
print_success "RGI detectado: ${RGI_VERSION}"
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
TOTAL_HITS=0

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
    SAMPLE_HITS=0
    
    for BIN in ${BINS}; do
        BIN_PROKKA="${SAMPLE_PROKKA}/${BIN}"
        BIN_OUT="${SAMPLE_OUT}/${BIN}"
        mkdir -p "${BIN_OUT}"
        
        # Buscar archivo de proteínas de Prokka
        PROTEIN_FILE="${BIN_PROKKA}/${BIN}.faa"
        
        if [ ! -f "${PROTEIN_FILE}" ]; then
            print_warning "${BIN}: No se encontró archivo de proteínas (.faa)"
            continue
        fi
        
        print_step "Analizando: ${BIN}"
        
        # Ejecutar RGI
        rgi main \
            --input_sequence "${PROTEIN_FILE}" \
            --output_file "${BIN_OUT}/${BIN}" \
            --input_type protein \
            --alignment_tool DIAMOND \
            --num_threads ${THREADS} \
            --clean \
            > "${BIN_OUT}/rgi.log" 2>&1
        
        if [ $? -eq 0 ]; then
            # Contar hits
            if [ -f "${BIN_OUT}/${BIN}.txt" ]; then
                HITS=$(tail -n +2 "${BIN_OUT}/${BIN}.txt" | wc -l)
                if [ ${HITS} -gt 0 ]; then
                    print_success "${HITS} genes de resistencia encontrados"
                    SAMPLE_HITS=$((SAMPLE_HITS + HITS))
                else
                    print_success "Sin genes de resistencia"
                fi
                SAMPLE_SUCCESS=$((SAMPLE_SUCCESS + 1))
                TOTAL_BINS=$((TOTAL_BINS + 1))
            else
                print_warning "RGI completado pero sin archivo de salida"
            fi
        else
            print_error "Error en RGI (ver log)"
        fi
    done
    
    if [ ${SAMPLE_SUCCESS} -gt 0 ]; then
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
        TOTAL_HITS=$((TOTAL_HITS + SAMPLE_HITS))
        echo ""
        print_success "${SAMPLE}: ${SAMPLE_SUCCESS} bins analizados, ${SAMPLE_HITS} genes de resistencia"
    else
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
    
    echo ""
done

# ==============================================================
# GENERAR RESUMEN COMBINADO
# ==============================================================

print_header "Generando Resumen Combinado"
echo ""

COMBINED="${OUTPUT_DIR}/RGI_All_Resistance_Genes.tsv"

{
    # Encabezado (tomar del primer archivo)
    FIRST_FILE=$(find "${OUTPUT_DIR}" -name "*.txt" | head -n1)
    if [ -n "${FIRST_FILE}" ]; then
        head -n1 "${FIRST_FILE}"
        # Datos de todos los archivos
        find "${OUTPUT_DIR}" -name "*.txt" -exec tail -n +2 {} \;
    fi
} > "${COMBINED}"

COMBINED_HITS=$(tail -n +2 "${COMBINED}" 2>/dev/null | wc -l)
print_success "${COMBINED_HITS} genes de resistencia totales"
echo ""

# ==============================================================
# RESUMEN FINAL
# ==============================================================

print_header "RGI COMPLETADO"
echo ""
echo "Resultados:"
echo -e "  ${GREEN}✓ Muestras exitosas: ${SUCCESS_COUNT}${NC}"
echo -e "  ${RED}✗ Muestras fallidas: ${FAIL_COUNT}${NC}"
echo -e "  ${CYAN}Total de bins analizados: ${TOTAL_BINS}${NC}"
echo -e "  ${YELLOW}Genes de resistencia: ${TOTAL_HITS}${NC}"
echo ""
echo "Resultados individuales en: ${OUTPUT_DIR}"
echo "Resumen combinado: ${COMBINED}"
echo ""
echo -e "${YELLOW}Siguiente paso:${NC}"
echo "  Ejecuta el módulo 9 (AntiSMASH) para analizar metabolitos secundarios:"
echo "  bash run_antismash.sh"
echo ""

exit 0
