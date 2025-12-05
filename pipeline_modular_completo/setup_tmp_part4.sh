#!/bin/bash

# ==============================================================
# Configuración de Directorios Temporales para /Part4
# ==============================================================
# 
# Este script configura los directorios temporales del pipeline
# para usar /Part4 en lugar de /tmp (root), evitando llenar el
# disco del sistema.
#
# Uso:
#   source setup_tmp_part4.sh
#   bash metagenomics_pipeline.sh
#
# ==============================================================

# Colores
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

echo "============================================================"
echo "  Configuración de Directorios Temporales"
echo "============================================================"
echo ""

# Directorio base para temporales
TMP_BASE="${TMP_BASE:-/mnt/Part4/Laboratory/tmp}"

# Verificar que el directorio base existe
if [ ! -d "/mnt/Part4/Laboratory" ]; then
    echo -e "${RED}✗ /mnt/Part4/Laboratory no encontrado${NC}"
    echo ""
    echo "Solución:"
    echo "  1. Verifica que /mnt/Part4 esté montado"
    echo "  2. Crea el directorio: mkdir -p /mnt/Part4/Laboratory/tmp"
    echo "  3. O cambia TMP_BASE a otra ubicación con espacio:"
    echo "     export TMP_BASE=/otra/ruta/tmp"
    echo "     source setup_tmp_part4.sh"
    echo ""
    return 1 2>/dev/null || exit 1
fi

# Crear estructura de directorios
echo "Creando directorios temporales en ${TMP_BASE}..."

DIRS=(
    "${TMP_BASE}/general"
    "${TMP_BASE}/megahit"
    "${TMP_BASE}/gtdbtk"
    "${TMP_BASE}/bowtie2"
    "${TMP_BASE}/concoct"
    "${TMP_BASE}/antismash"
    "${TMP_BASE}/maxbin2"
    "${TMP_BASE}/prokka"
)

for dir in "${DIRS[@]}"; do
    if mkdir -p "${dir}" 2>/dev/null; then
        echo -e "  ${GREEN}✓${NC} ${dir}"
    else
        echo -e "  ${RED}✗${NC} ${dir} (sin permisos)"
    fi
done

echo ""

# Configurar permisos (si es posible)
chmod -R 1777 "${TMP_BASE}" 2>/dev/null || true

# Exportar variables de entorno
export TMPDIR="${TMP_BASE}/general"
export TEMP="${TMPDIR}"
export TMP="${TMPDIR}"

# Variables específicas por herramienta
export MEGAHIT_TMP="${TMP_BASE}/megahit"
export GTDBTK_TMP="${TMP_BASE}/gtdbtk"
export BOWTIE2_TMP="${TMP_BASE}/bowtie2"
export CONCOCT_TMP="${TMP_BASE}/concoct"
export ANTISMASH_TMP="${TMP_BASE}/antismash"
export MAXBIN2_TMP="${TMP_BASE}/maxbin2"
export PROKKA_TMP="${TMP_BASE}/prokka"

# Verificar espacio disponible
FREE_SPACE=$(df -h --output=avail "${TMP_BASE}" | tail -n1 | tr -d ' ')
FREE_SPACE_GB=$(df -BG --output=avail "${TMP_BASE}" | tail -n1 | tr -dc '0-9')

echo "Configuración:"
echo "  Directorio base: ${TMP_BASE}"
echo "  Espacio disponible: ${FREE_SPACE}"
echo ""

# Advertir si hay poco espacio
if [ "${FREE_SPACE_GB}" -lt 500 ]; then
    echo -e "${YELLOW}⚠ Advertencia: Solo ${FREE_SPACE} disponibles${NC}"
    echo "  Se recomienda >500 GB para el pipeline completo"
    echo ""
fi

# Mostrar variables exportadas
echo "Variables de entorno configuradas:"
echo "  TMPDIR=${TMPDIR}"
echo "  MEGAHIT_TMP=${MEGAHIT_TMP}"
echo "  GTDBTK_TMP=${GTDBTK_TMP}"
echo "  CONCOCT_TMP=${CONCOCT_TMP}"
echo "  ANTISMASH_TMP=${ANTISMASH_TMP}"
echo ""

echo -e "${GREEN}✓ Configuración completada${NC}"
echo ""
echo "Siguiente paso:"
echo "  bash metagenomics_pipeline.sh"
echo ""
echo "Para limpiar archivos temporales después:"
echo "  rm -rf ${TMP_BASE}/*"
echo ""
