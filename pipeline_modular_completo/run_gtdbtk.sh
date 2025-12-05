#!/bin/bash

# ==============================================================
# Clasificación taxonómica de bins con GTDB-Tk v2.5.2
# Base de datos: r226 (con soporte Skani)
# ==============================================================

# Detectar directorio base (donde está el script)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"

# Directorios (usa variables de entorno o valores por defecto)
BINNING_DIR="${BINNING_DIR:-${PROJECT_DIR}/output/binning}"
OUTPUT_DIR="${OUTPUT_DIR:-${PROJECT_DIR}/output/gtdbtk}"
THREADS="${THREADS:-40}"

# Selección de bins (usa variable de entorno o auto)
# Opciones: dastool, metabat2, maxbin2, concoct, auto
BINS_SOURCE="${BINS_SOURCE:-auto}"

# Configuración de GTDB-Tk
# Intentar detectar automáticamente si no está configurado
if [ -z "${GTDBTK_DATA_PATH}" ]; then
    # Buscar en ubicaciones comunes
    POSSIBLE_PATHS=(
        "/data/database/gtdbtk_251103/gtdbtk_data_release226"
        "/data/database/gtdbtk_data_release226"
        "/files/database/gtdbtk_data_release226"
        "${HOME}/database/gtdbtk_data_release226"
        "/opt/database/gtdbtk_data_release226"
    )
    
    for path in "${POSSIBLE_PATHS[@]}"; do
        if [ -d "${path}" ]; then
            export GTDBTK_DATA_PATH="${path}"
            break
        fi
    done
fi

# Si aún no está configurado, usar valor por defecto
export GTDBTK_DATA_PATH="${GTDBTK_DATA_PATH:-/data/database/gtdbtk_251103/gtdbtk_data_release226}"

# Configuración de Skani (NUEVO en v2.5+)
# Skani es mucho más rápido que ANI para clasificación
# Si tu base de datos tiene la carpeta skani/, déjalo habilitado
if [ -d "${GTDBTK_DATA_PATH}/skani" ]; then
    export GTDBTK_DISABLE_SKANI=0
    SKANI_STATUS="✓ Habilitado (rápido)"
else
    export GTDBTK_DISABLE_SKANI=1
    SKANI_STATUS="⚠ Deshabilitado (base de datos sin skani/)"
fi

# Directorio temporal
export TMPDIR="${TMPDIR:-/tmp/gtdbtk_$$}"

# Colores
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
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

print_header "GTDB-Tk - Clasificación Taxonómica"
echo ""
echo "Configuración:"
echo "  Base de datos: ${GTDBTK_DATA_PATH}"
echo "  Directorio de bins: ${BINNING_DIR}"
echo "  Directorio de salida: ${OUTPUT_DIR}"
echo "  Fuente de bins: ${BINS_SOURCE}"
echo "  Threads: ${THREADS}"
echo "  Skani: ${SKANI_STATUS}"
echo "  TMPDIR: ${TMPDIR}"
echo ""

# Crear directorios
mkdir -p "${TMPDIR}"
mkdir -p "${OUTPUT_DIR}"

# Verificar espacio en disco
FREE_SPACE=$(df -h --output=avail -BG "${TMPDIR}" | tail -n1 | tr -dc '0-9')
if [ "${FREE_SPACE}" -lt 150 ]; then
    print_warning "Solo ${FREE_SPACE} GB libres en ${TMPDIR}. Se recomienda >150 GB."
else
    print_success "Espacio temporal disponible: ${FREE_SPACE} GB"
fi
echo ""

# Verificar que GTDB-Tk esté disponible
if ! command -v gtdbtk &> /dev/null; then
    print_error "GTDB-Tk no está disponible"
    echo ""
    echo "Solución:"
    echo "  1. Activa el ambiente: micromamba activate gtdbtk"
    echo "  2. Verifica la instalación: gtdbtk --version"
    echo ""
    exit 1
fi

# Verificar versión de GTDB-Tk
GTDBTK_VERSION=$(gtdbtk --version 2>&1 | grep -oP 'v\K[0-9.]+' | head -n1)
print_success "GTDB-Tk v${GTDBTK_VERSION} detectado"

# Verificar base de datos
if [ -z "${GTDBTK_DATA_PATH}" ]; then
    print_error "Variable GTDBTK_DATA_PATH no configurada"
    echo ""
    echo "Solución:"
    echo "  1. Encuentra tu base de datos GTDB-Tk:"
    echo "     find /data /files /opt -name 'gtdbtk_data_release*' -type d 2>/dev/null"
    echo ""
    echo "  2. Configura la variable:"
    echo "     export GTDBTK_DATA_PATH=/ruta/a/gtdbtk_data_release226"
    echo ""
    echo "  3. O descarga la base de datos r226:"
    echo "     wget https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_r226_data.tar.gz"
    echo "     tar -xzf gtdbtk_r226_data.tar.gz"
    echo ""
    exit 1
fi

if [ ! -d "${GTDBTK_DATA_PATH}" ]; then
    print_error "Base de datos no encontrada: ${GTDBTK_DATA_PATH}"
    echo ""
    echo "El directorio no existe. Verifica la ruta o descarga la base de datos."
    echo ""
    exit 1
fi

DB_VERSION=$(basename "${GTDBTK_DATA_PATH}")
print_success "Base de datos: ${DB_VERSION}"
echo ""

# ==============================================================
# PROCESAMIENTO
# ==============================================================

# Encontrar muestras
SAMPLES=$(ls -d ${BINNING_DIR}/*/ 2>/dev/null | xargs -n 1 basename)

if [ -z "${SAMPLES}" ]; then
    print_error "No se encontraron muestras en ${BINNING_DIR}"
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

for SAMPLE in ${SAMPLES}; do
    print_header "Procesando: ${SAMPLE}"
    
    SAMPLE_DIR="${BINNING_DIR}/${SAMPLE}"
    SAMPLE_OUT="${OUTPUT_DIR}/${SAMPLE}"
    
    # Seleccionar bins según configuración
    BINS_DIR=""
    BINS_EXT=""
    BINS_NAME=""
    
    case "${BINS_SOURCE}" in
        dastool)
            if [ -d "${SAMPLE_DIR}/dastool/DASTool_DASTool_bins" ]; then
                BINS_DIR="${SAMPLE_DIR}/dastool/DASTool_DASTool_bins"
                BINS_EXT="fa"
                BINS_NAME="DAS Tool (refinados)"
            fi
            ;;
        metabat2)
            if [ -d "${SAMPLE_DIR}/metabat2" ]; then
                BINS_DIR="${SAMPLE_DIR}/metabat2"
                BINS_EXT="fa"
                BINS_NAME="MetaBAT2"
            fi
            ;;
        maxbin2)
            if [ -d "${SAMPLE_DIR}/maxbin2" ]; then
                BINS_DIR="${SAMPLE_DIR}/maxbin2"
                BINS_EXT="fasta"
                BINS_NAME="MaxBin2"
            fi
            ;;
        concoct)
            if [ -d "${SAMPLE_DIR}/concoct/bins" ]; then
                BINS_DIR="${SAMPLE_DIR}/concoct/bins"
                BINS_EXT="fa"
                BINS_NAME="CONCOCT"
            fi
            ;;
        auto|*)
            # Buscar en orden de preferencia
            if [ -d "${SAMPLE_DIR}/dastool/DASTool_DASTool_bins" ]; then
                BINS_DIR="${SAMPLE_DIR}/dastool/DASTool_DASTool_bins"
                BINS_EXT="fa"
                BINS_NAME="DAS Tool (refinados)"
            elif [ -d "${SAMPLE_DIR}/metabat2" ]; then
                BINS_DIR="${SAMPLE_DIR}/metabat2"
                BINS_EXT="fa"
                BINS_NAME="MetaBAT2"
            elif [ -d "${SAMPLE_DIR}/maxbin2" ]; then
                BINS_DIR="${SAMPLE_DIR}/maxbin2"
                BINS_EXT="fasta"
                BINS_NAME="MaxBin2"
            elif [ -d "${SAMPLE_DIR}/concoct/bins" ]; then
                BINS_DIR="${SAMPLE_DIR}/concoct/bins"
                BINS_EXT="fa"
                BINS_NAME="CONCOCT"
            fi
            ;;
    esac
    
    # Verificar que se encontraron bins
    if [ -z "${BINS_DIR}" ]; then
        if [ "${BINS_SOURCE}" = "auto" ]; then
            print_error "No se encontraron bins para ${SAMPLE}"
        else
            print_error "No se encontraron bins de ${BINS_SOURCE} para ${SAMPLE}"
        fi
        FAIL_COUNT=$((FAIL_COUNT + 1))
        echo ""
        continue
    fi
    
    # Contar bins
    BIN_COUNT=$(ls ${BINS_DIR}/*.${BINS_EXT} 2>/dev/null | wc -l)
    
    if [ ${BIN_COUNT} -eq 0 ]; then
        print_error "No hay bins (*.${BINS_EXT}) en ${BINS_DIR}"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        echo ""
        continue
    fi
    
    print_step "Fuente de bins: ${BINS_NAME}"
    print_step "Bins encontrados: ${BIN_COUNT}"
    print_step "Directorio: ${BINS_DIR}"
    echo ""
    
    # Crear directorio de salida
    mkdir -p "${SAMPLE_OUT}"
    
    # Verificar si ya fue procesado
    if [ -f "${SAMPLE_OUT}/gtdbtk.bac120.summary.tsv" ] || [ -f "${SAMPLE_OUT}/gtdbtk.ar53.summary.tsv" ]; then
        print_warning "Resultados previos detectados, omitiendo"
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
        echo ""
        continue
    fi
    
    # Ejecutar GTDB-Tk
    print_step "Ejecutando GTDB-Tk classify_wf..."
    echo ""
    
    # Opciones según disponibilidad de Skani
    if [ "${GTDBTK_DISABLE_SKANI}" -eq 0 ]; then
        # Con Skani (más rápido)
        gtdbtk classify_wf \
            --genome_dir "${BINS_DIR}" \
            --out_dir "${SAMPLE_OUT}" \
            --extension "${BINS_EXT}" \
            --cpus ${THREADS} \
            --skip_ani_screen \
            > "${SAMPLE_OUT}/gtdbtk.log" 2>&1
    else
        # Sin Skani (compatible con bases antiguas)
        gtdbtk classify_wf \
            --genome_dir "${BINS_DIR}" \
            --out_dir "${SAMPLE_OUT}" \
            --extension "${BINS_EXT}" \
            --cpus ${THREADS} \
            --skip_ani_screen \
            > "${SAMPLE_OUT}/gtdbtk.log" 2>&1
    fi
    
    STATUS=$?
    
    if [ ${STATUS} -eq 0 ]; then
        print_success "Clasificación completada"
        
        # Contar bins clasificados
        BAC_COUNT=0
        AR_COUNT=0
        
        if [ -f "${SAMPLE_OUT}/gtdbtk.bac120.summary.tsv" ]; then
            BAC_COUNT=$(tail -n +2 "${SAMPLE_OUT}/gtdbtk.bac120.summary.tsv" | wc -l)
        fi
        
        if [ -f "${SAMPLE_OUT}/gtdbtk.ar53.summary.tsv" ]; then
            AR_COUNT=$(tail -n +2 "${SAMPLE_OUT}/gtdbtk.ar53.summary.tsv" | wc -l)
        fi
        
        echo ""
        print_step "Resultados:"
        echo -e "    ${GREEN}Bacterias:${NC} ${BAC_COUNT} bins"
        echo -e "    ${GREEN}Arqueas:${NC} ${AR_COUNT} bins"
        
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    else
        print_error "Error en clasificación (ver log: ${SAMPLE_OUT}/gtdbtk.log)"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
    
    echo ""
done

# ==============================================================
# COMBINAR RESULTADOS
# ==============================================================

print_header "Generando Resumen Combinado"
echo ""

COMBINED_BAC="${OUTPUT_DIR}/GTDBTK_All_Bacteria.tsv"
COMBINED_AR="${OUTPUT_DIR}/GTDBTK_All_Archaea.tsv"

# Combinar bacterias
print_step "Combinando clasificaciones bacterianas..."
{
    # Encabezado (tomar del primer archivo)
    FIRST_BAC=$(find "${OUTPUT_DIR}" -name "gtdbtk.bac120.summary.tsv" | head -n1)
    if [ -n "${FIRST_BAC}" ]; then
        head -n1 "${FIRST_BAC}"
        # Datos de todos los archivos
        find "${OUTPUT_DIR}" -name "gtdbtk.bac120.summary.tsv" -exec tail -n +2 {} \;
    fi
} > "${COMBINED_BAC}"

BAC_TOTAL=$(tail -n +2 "${COMBINED_BAC}" 2>/dev/null | wc -l)
print_success "${BAC_TOTAL} bins bacterianos clasificados"

# Combinar arqueas
print_step "Combinando clasificaciones arqueales..."
{
    # Encabezado (tomar del primer archivo)
    FIRST_AR=$(find "${OUTPUT_DIR}" -name "gtdbtk.ar53.summary.tsv" | head -n1)
    if [ -n "${FIRST_AR}" ]; then
        head -n1 "${FIRST_AR}"
        # Datos de todos los archivos
        find "${OUTPUT_DIR}" -name "gtdbtk.ar53.summary.tsv" -exec tail -n +2 {} \;
    fi
} > "${COMBINED_AR}"

AR_TOTAL=$(tail -n +2 "${COMBINED_AR}" 2>/dev/null | wc -l)
print_success "${AR_TOTAL} bins arqueales clasificados"

echo ""

# ==============================================================
# RESUMEN FINAL
# ==============================================================

print_header "GTDB-Tk COMPLETADO"
echo ""
echo "Resultados:"
echo -e "  ${GREEN}✓ Exitosas: ${SUCCESS_COUNT}${NC}"
echo -e "  ${RED}✗ Fallidas: ${FAIL_COUNT}${NC}"
echo ""
echo "Clasificaciones:"
echo -e "  ${CYAN}Bacterias:${NC} ${BAC_TOTAL} bins → ${COMBINED_BAC}"
echo -e "  ${CYAN}Arqueas:${NC} ${AR_TOTAL} bins → ${COMBINED_AR}"
echo ""
echo "Archivos individuales en: ${OUTPUT_DIR}/"
echo ""

# Limpiar directorio temporal
rm -rf "${TMPDIR}"

exit 0
