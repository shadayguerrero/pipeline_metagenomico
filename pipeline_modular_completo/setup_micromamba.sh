#!/bin/bash

# Script para configurar micromamba en el PATH
# Ejecuta este script antes de usar el pipeline:
# source setup_micromamba.sh

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  Configuración de Micromamba${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""

# Buscar micromamba en ubicaciones comunes
MICROMAMBA_LOCATIONS=(
    "${HOME}/micromamba/bin/micromamba"
    "${HOME}/.local/bin/micromamba"
    "/home_local/camda/micromamba/bin/micromamba"
    "/usr/local/bin/micromamba"
    "/opt/micromamba/bin/micromamba"
)

MICROMAMBA_FOUND=""

for location in "${MICROMAMBA_LOCATIONS[@]}"; do
    if [ -f "${location}" ]; then
        MICROMAMBA_FOUND="${location}"
        break
    fi
done

if [ -z "${MICROMAMBA_FOUND}" ]; then
    echo -e "${RED}✗ No se encontró micromamba en las ubicaciones comunes${NC}"
    echo ""
    echo "Ubicaciones buscadas:"
    for location in "${MICROMAMBA_LOCATIONS[@]}"; do
        echo "  - ${location}"
    done
    echo ""
    echo -e "${YELLOW}Por favor, instala micromamba o especifica su ubicación:${NC}"
    echo ""
    echo "  # Para instalar micromamba:"
    echo "  \"\${SHELL}\" <(curl -L micro.mamba.pm/install.sh)"
    echo ""
    echo "  # O añade manualmente al PATH:"
    echo "  export PATH=\"/ruta/a/micromamba/bin:\$PATH\""
    echo ""
    exit 1
fi

# Obtener el directorio bin
MICROMAMBA_BIN_DIR=$(dirname "${MICROMAMBA_FOUND}")

echo -e "${GREEN}✓ Micromamba encontrado en: ${MICROMAMBA_FOUND}${NC}"
echo ""

# Añadir al PATH si no está ya
if [[ ":$PATH:" != *":${MICROMAMBA_BIN_DIR}:"* ]]; then
    export PATH="${MICROMAMBA_BIN_DIR}:$PATH"
    echo -e "${GREEN}✓ Añadido al PATH: ${MICROMAMBA_BIN_DIR}${NC}"
else
    echo -e "${GREEN}✓ Ya está en el PATH${NC}"
fi

# Inicializar el shell hook
eval "$(micromamba shell hook --shell bash)"

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ Shell hook inicializado${NC}"
else
    echo -e "${RED}✗ Error al inicializar shell hook${NC}"
    exit 1
fi

echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  Configuración Completada${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "Ahora puedes ejecutar el pipeline:"
echo "  bash metagenomics_pipeline.sh"
echo ""
echo "O activar ambientes manualmente:"
echo "  micromamba activate kraken2"
echo ""

# Listar ambientes disponibles
echo "Ambientes disponibles:"
micromamba env list 2>/dev/null | grep -v "^#" | grep -v "^$"
echo ""

