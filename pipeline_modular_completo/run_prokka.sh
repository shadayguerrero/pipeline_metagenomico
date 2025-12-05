#!/bin/bash

# ==============================================================
# Prokka: anotación funcional de bins según clasificación GTDB-Tk
# ==============================================================

# Detectar directorio base (donde está el script)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"

# Directorios (usa variables de entorno o valores por defecto)
BINNING_DIR="${BINNING_DIR:-${PROJECT_DIR}/output/binning}"
GTDBTK_DIR="${GTDBTK_DIR:-${PROJECT_DIR}/output/gtdbtk}"
OUTPUT_DIR="${OUTPUT_DIR:-${PROJECT_DIR}/output/prokka}"
THREADS="${THREADS:-8}"

# Selección de bins (usa variable de entorno o auto)
# Opciones: dastool, metabat2, maxbin2, concoct, auto
BINS_SOURCE="${BINS_SOURCE:-auto}"

# Parámetros de Prokka
MIN_CONTIG_LEN=200
USE_GENUS=0   # Cambia a 1 si quieres usar --usegenus --genus/--species

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

# Función de anotación por bin
annotate_bin() {
    local bin_fa="$1"
    local sample="$2"
    local bin_name="$3"
    local kingdom="$4"
    local genus="$5"
    local species="$6"
    
    local out_bin="${OUTPUT_DIR}/${sample}/${bin_name}"
    mkdir -p "${out_bin}"
    
    print_step "Anotando: ${bin_name} (${kingdom})"
    
    extra_flags=()
    if [[ "${USE_GENUS}" -eq 1 && "${genus}" != "-" ]]; then
        extra_flags+=( --usegenus --genus "${genus}" )
        [[ "${species}" != "-" ]] && extra_flags+=( --species "${species}" )
    fi
    
    prokka \
        --outdir "${out_bin}" \
        --prefix "${bin_name}" \
        --cpus "${THREADS}" \
        --mincontiglen "${MIN_CONTIG_LEN}" \
        --metagenome \
        --force \
        --kingdom "${kingdom}" \
        "${extra_flags[@]}" \
        "${bin_fa}" > "${out_bin}/prokka.log" 2>&1
    
    if [ $? -eq 0 ]; then
        print_success "Anotación completada"
    else
        print_error "Error en anotación (ver log)"
    fi
}

# ==============================================================
# PREPARACIÓN
# ==============================================================

print_header "Prokka - Anotación Funcional"
echo ""
echo "Configuración:"
echo "  Directorio de bins: ${BINNING_DIR}"
echo "  Directorio GTDB-Tk: ${GTDBTK_DIR}"
echo "  Directorio de salida: ${OUTPUT_DIR}"
echo "  Fuente de bins: ${BINS_SOURCE}"
echo "  Threads: ${THREADS}"
echo "  Longitud mínima de contigs: ${MIN_CONTIG_LEN} bp"
echo ""

# Verificar que Prokka esté disponible
if ! command -v prokka &> /dev/null; then
    print_error "Prokka no está disponible"
    echo ""
    echo "Solución:"
    echo "  1. Activa el ambiente: micromamba activate prokka"
    echo "  2. Verifica la instalación: prokka --version"
    echo ""
    exit 1
fi

PROKKA_VERSION=$(prokka --version 2>&1 | head -n1)
print_success "Prokka detectado: ${PROKKA_VERSION}"
echo ""

# Crear directorio de salida
mkdir -p "${OUTPUT_DIR}"

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

# ==============================================================
# PROCESAMIENTO
# ==============================================================

SUCCESS_COUNT=0
FAIL_COUNT=0
TOTAL_BINS=0

for SAMPLE in ${SAMPLES}; do
    print_header "Procesando: ${SAMPLE}"
    
    SAMPLE_DIR="${BINNING_DIR}/${SAMPLE}"
    GTDBTK_SAMPLE="${GTDBTK_DIR}/${SAMPLE}"
    
    # Seleccionar bins según configuración
    BINS_DIR=""
    BINS_EXT=""
    BINS_NAME=""
    
    case "${BINS_SOURCE}" in
        dastool)
            if [ -d "${SAMPLE_DIR}/dastool/DASTool_DASTool_bins" ]; then
                BINS_DIR="${SAMPLE_DIR}/dastool/DASTool_DASTool_bins"
                BINS_EXT="fa"
                BINS_NAME="DAS Tool"
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
                BINS_NAME="DAS Tool"
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
    
    # Verificar bins
    if [ -z "${BINS_DIR}" ]; then
        print_error "No se encontraron bins de ${BINS_SOURCE}"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        echo ""
        continue
    fi
    
    BIN_COUNT=$(ls ${BINS_DIR}/*.${BINS_EXT} 2>/dev/null | wc -l)
    
    if [ ${BIN_COUNT} -eq 0 ]; then
        print_error "No hay bins en ${BINS_DIR}"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        echo ""
        continue
    fi
    
    print_step "Fuente de bins: ${BINS_NAME}"
    print_step "Bins encontrados: ${BIN_COUNT}"
    echo ""
    
    # Verificar resultados de GTDB-Tk
    SUM_BAC="${GTDBTK_SAMPLE}/gtdbtk.bac120.summary.tsv"
    SUM_ARC="${GTDBTK_SAMPLE}/gtdbtk.ar53.summary.tsv"
    
    if [ ! -f "${SUM_BAC}" ] && [ ! -f "${SUM_ARC}" ]; then
        print_warning "No se encontraron resultados de GTDB-Tk para ${SAMPLE}"
        print_warning "Ejecuta primero: bash run_gtdbtk.sh"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        echo ""
        continue
    fi
    
    # Crear mapa de clasificaciones
    TMP_MAP=$(mktemp)
    
    # Bacterias
    if [ -f "${SUM_BAC}" ]; then
        awk -F'\t' 'NR>1{
            cmd="basename \""$1"\""; cmd|getline b; close(cmd);
            sub(/\.[Ff][Aa][Ss]?[Tt]?[Aa]?$/, "", b);
            cl=$2; g="-"; s="-";
            n=split(cl,a,";");
            for(i=1;i<=n;i++){ 
                if(a[i]~/^g__/){g=a[i];sub(/^g__/,"",g)}; 
                if(a[i]~/^s__/){s=a[i];sub(/^s__/,"",s)} 
            }
            print b "\tBacteria\t" g "\t" s
        }' "${SUM_BAC}" >> "${TMP_MAP}"
    fi
    
    # Arqueas
    if [ -f "${SUM_ARC}" ]; then
        awk -F'\t' 'NR>1{
            cmd="basename \""$1"\""; cmd|getline b; close(cmd);
            sub(/\.[Ff][Aa][Ss]?[Tt]?[Aa]?$/, "", b);
            cl=$2; g="-"; s="-";
            n=split(cl,a,";");
            for(i=1;i<=n;i++){ 
                if(a[i]~/^g__/){g=a[i];sub(/^g__/,"",g)}; 
                if(a[i]~/^s__/){s=a[i];sub(/^s__/,"",s)} 
            }
            print b "\tArchaea\t" g "\t" s
        }' "${SUM_ARC}" >> "${TMP_MAP}"
    fi
    
    # Procesar cada bin
    SAMPLE_SUCCESS=0
    
    for bin_fa in "${BINS_DIR}"/*.${BINS_EXT}; do
        [ -f "${bin_fa}" ] || continue
        
        bin_base=$(basename "${bin_fa}")
        bin_name="${bin_base%.*}"
        
        # Buscar clasificación
        line=$(grep -F -m1 -w "${bin_name}" "${TMP_MAP}" || true)
        
        if [ -z "${line}" ]; then
            print_warning "${bin_name} no tiene clasificación GTDB-Tk, omitiendo"
            continue
        fi
        
        kingdom=$(echo "${line}" | awk -F'\t' '{print $2}')
        genus=$(echo "${line}" | awk -F'\t' '{print $3}')
        species=$(echo "${line}" | awk -F'\t' '{print $4}')
        
        if [ "${kingdom}" != "Bacteria" ] && [ "${kingdom}" != "Archaea" ]; then
            print_warning "${bin_name}: reino ${kingdom} no soportado, omitiendo"
            continue
        fi
        
        annotate_bin "${bin_fa}" "${SAMPLE}" "${bin_name}" "${kingdom}" "${genus}" "${species}"
        SAMPLE_SUCCESS=$((SAMPLE_SUCCESS + 1))
        TOTAL_BINS=$((TOTAL_BINS + 1))
    done
    
    rm -f "${TMP_MAP}"
    
    if [ ${SAMPLE_SUCCESS} -gt 0 ]; then
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
        print_success "${SAMPLE_SUCCESS} bins anotados"
    else
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
    
    echo ""
done

# ==============================================================
# RESUMEN FINAL
# ==============================================================

print_header "PROKKA COMPLETADO"
echo ""
echo "Resultados:"
echo -e "  ${GREEN}✓ Muestras exitosas: ${SUCCESS_COUNT}${NC}"
echo -e "  ${RED}✗ Muestras fallidas: ${FAIL_COUNT}${NC}"
echo -e "  ${CYAN}Total de bins anotados: ${TOTAL_BINS}${NC}"
echo ""
echo "Anotaciones en: ${OUTPUT_DIR}"
echo ""
echo -e "${YELLOW}Siguiente paso:${NC}"
echo "  Ejecuta el módulo 8 (RGI) para analizar resistencia a antibióticos:"
echo "  bash run_rgi.sh"
echo ""

exit 0
