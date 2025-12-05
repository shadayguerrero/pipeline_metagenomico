#!/bin/bash

# ============================================================================
# Script de Instalación de Ambientes para Módulos 6-9
# Pipeline Metagenómico Modular
# ============================================================================

set -e

echo "============================================================================"
echo "  Instalación de Ambientes - Módulos 6-9"
echo "============================================================================"
echo ""

# Colores
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

# Verificar que micromamba esté disponible
if ! command -v micromamba &> /dev/null; then
    echo -e "${RED}Error: micromamba no está disponible${NC}"
    echo "Instala micromamba desde: https://mamba.readthedocs.io/en/latest/installation.html"
    exit 1
fi

echo -e "${GREEN}✓ micromamba encontrado${NC}"
echo ""

# ============================================================================
# MÓDULO 6: GTDB-Tk
# ============================================================================

echo -e "${YELLOW}[1/4] Instalando ambiente para GTDB-Tk...${NC}"

if micromamba env list | grep -q "^gtdbtk "; then
    echo -e "${YELLOW}⚠ Ambiente 'gtdbtk' ya existe. ¿Deseas reinstalarlo? (s/n)${NC}"
    read -r response
    if [[ "$response" =~ ^[Ss]$ ]]; then
        micromamba env remove -n gtdbtk -y
        micromamba create -n gtdbtk -c bioconda -c conda-forge gtdbtk -y
    else
        echo -e "${GREEN}✓ Usando ambiente existente${NC}"
    fi
else
    micromamba create -n gtdbtk -c bioconda -c conda-forge gtdbtk -y
fi

echo -e "${GREEN}✓ Ambiente 'gtdbtk' instalado${NC}"
echo ""

# Verificar base de datos GTDB
echo -e "${YELLOW}Verificando base de datos GTDB...${NC}"
if [ -z "${GTDBTK_DATA_PATH}" ]; then
    echo -e "${RED}⚠ Variable GTDBTK_DATA_PATH no configurada${NC}"
    echo "Descarga la base de datos desde: https://data.gtdb.ecogenomic.org/"
    echo ""
    echo "Comandos sugeridos:"
    echo "  wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz"
    echo "  tar -xzf gtdbtk_data.tar.gz"
    echo "  export GTDBTK_DATA_PATH=/path/to/gtdbtk_data"
    echo ""
else
    if [ -d "${GTDBTK_DATA_PATH}" ]; then
        echo -e "${GREEN}✓ Base de datos GTDB encontrada en: ${GTDBTK_DATA_PATH}${NC}"
    else
        echo -e "${RED}✗ Directorio no encontrado: ${GTDBTK_DATA_PATH}${NC}"
    fi
fi
echo ""

# ============================================================================
# MÓDULO 7: Prokka
# ============================================================================

echo -e "${YELLOW}[2/4] Instalando ambiente para Prokka...${NC}"

if micromamba env list | grep -q "^prokka "; then
    echo -e "${YELLOW}⚠ Ambiente 'prokka' ya existe. ¿Deseas reinstalarlo? (s/n)${NC}"
    read -r response
    if [[ "$response" =~ ^[Ss]$ ]]; then
        micromamba env remove -n prokka -y
        micromamba create -n prokka -c bioconda -c conda-forge prokka -y
    else
        echo -e "${GREEN}✓ Usando ambiente existente${NC}"
    fi
else
    micromamba create -n prokka -c bioconda -c conda-forge prokka -y
fi

echo -e "${GREEN}✓ Ambiente 'prokka' instalado${NC}"
echo ""

# Actualizar bases de datos de Prokka
echo -e "${YELLOW}Actualizando bases de datos de Prokka...${NC}"
eval "$(micromamba shell hook --shell bash)"
micromamba activate prokka
prokka --setupdb || echo -e "${YELLOW}⚠ No se pudo actualizar la base de datos (puede ser normal)${NC}"
micromamba deactivate
echo ""

# ============================================================================
# MÓDULO 8: RGI
# ============================================================================

echo -e "${YELLOW}[3/4] Instalando ambiente para RGI...${NC}"

if micromamba env list | grep -q "^rgi "; then
    echo -e "${YELLOW}⚠ Ambiente 'rgi' ya existe. ¿Deseas reinstalarlo? (s/n)${NC}"
    read -r response
    if [[ "$response" =~ ^[Ss]$ ]]; then
        micromamba env remove -n rgi -y
        micromamba create -n rgi -c bioconda -c conda-forge rgi -y
    else
        echo -e "${GREEN}✓ Usando ambiente existente${NC}"
    fi
else
    micromamba create -n rgi -c bioconda -c conda-forge rgi -y
fi

echo -e "${GREEN}✓ Ambiente 'rgi' instalado${NC}"
echo ""

# Descargar base de datos CARD
echo -e "${YELLOW}Descargando base de datos CARD...${NC}"
eval "$(micromamba shell hook --shell bash)"
micromamba activate rgi

# Verificar si CARD ya está cargada
if rgi database --version 2>/dev/null | grep -q "CARD"; then
    echo -e "${GREEN}✓ Base de datos CARD ya está cargada${NC}"
else
    echo "Descargando CARD..."
    wget -O card_data.tar.bz2 https://card.mcmaster.ca/latest/data
    tar -xjf card_data.tar.bz2
    rgi load --card_json card.json
    rm -f card_data.tar.bz2
    echo -e "${GREEN}✓ Base de datos CARD descargada y cargada${NC}"
fi

micromamba deactivate
echo ""

# ============================================================================
# MÓDULO 9: AntiSMASH
# ============================================================================

echo -e "${YELLOW}[4/4] Instalando ambiente para AntiSMASH...${NC}"

if micromamba env list | grep -q "^antismash "; then
    echo -e "${YELLOW}⚠ Ambiente 'antismash' ya existe. ¿Deseas reinstalarlo? (s/n)${NC}"
    read -r response
    if [[ "$response" =~ ^[Ss]$ ]]; then
        micromamba env remove -n antismash -y
        micromamba create -n antismash -c bioconda -c conda-forge antismash -y
    else
        echo -e "${GREEN}✓ Usando ambiente existente${NC}"
    fi
else
    micromamba create -n antismash -c bioconda -c conda-forge antismash -y
fi

echo -e "${GREEN}✓ Ambiente 'antismash' instalado${NC}"
echo ""

# Descargar bases de datos de AntiSMASH
echo -e "${YELLOW}Descargando bases de datos de AntiSMASH...${NC}"
echo "Nota: Las bases de datos se descargarán automáticamente en la primera ejecución"
echo "Ubicación: ~/.local/share/antismash/"
echo ""

# ============================================================================
# RESUMEN
# ============================================================================

echo ""
echo "============================================================================"
echo -e "${GREEN}  ✓ Instalación Completada${NC}"
echo "============================================================================"
echo ""
echo "Ambientes instalados:"
echo "  1. gtdbtk   - Taxonomía de bins"
echo "  2. prokka   - Anotación funcional"
echo "  3. rgi      - Resistencia a antibióticos"
echo "  4. antismash - Metabolitos secundarios"
echo ""
echo "Verificar instalación:"
echo "  micromamba env list"
echo ""
echo "Activar un ambiente:"
echo "  micromamba activate gtdbtk"
echo ""
echo "Bases de datos requeridas:"
echo "  ✓ CARD (RGI) - Instalada"
echo "  ⚠ GTDB (GTDB-Tk) - Verifica GTDBTK_DATA_PATH"
echo "  ⚠ AntiSMASH - Se descarga automáticamente"
echo ""
echo "============================================================================"
