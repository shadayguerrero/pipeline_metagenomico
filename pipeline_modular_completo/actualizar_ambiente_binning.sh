#!/bin/bash

# Script para actualizar el ambiente de binning con todas las herramientas
# Versión: 2.0

# Colores
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Actualización del Ambiente de Binning${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Verificar que micromamba esté disponible
if ! command -v micromamba &> /dev/null; then
    echo -e "${RED}✗ Error: micromamba no está disponible${NC}"
    echo ""
    echo "Solución:"
    echo "  1. Instala micromamba desde: https://mamba.readthedocs.io/en/latest/installation.html"
    echo "  2. O agrega al PATH: export PATH=\"/ruta/a/micromamba/bin:\$PATH\""
    echo ""
    exit 1
fi

echo -e "${GREEN}✓ micromamba encontrado${NC}"
echo ""

# Verificar si el ambiente existe
if micromamba env list | grep -q "^binning "; then
    echo -e "${YELLOW}⚠ El ambiente 'binning' ya existe${NC}"
    echo ""
    read -p "¿Deseas reinstalarlo? (s/n): " -n 1 -r
    echo ""
    
    if [[ $REPLY =~ ^[Ss]$ ]]; then
        echo -e "${CYAN}Eliminando ambiente anterior...${NC}"
        micromamba env remove -n binning -y
        echo -e "${GREEN}✓ Ambiente eliminado${NC}"
        echo ""
    else
        echo -e "${CYAN}Actualizando ambiente existente...${NC}"
        echo ""
        
        # Activar ambiente
        eval "$(micromamba shell hook --shell bash)"
        micromamba activate binning
        
        # Instalar herramientas faltantes
        echo -e "${CYAN}Instalando herramientas adicionales...${NC}"
        micromamba install -y -c bioconda -c conda-forge \
            maxbin2 \
            concoct \
            das_tool \
            bbmap \
            diamond
        
        if [ $? -eq 0 ]; then
            echo ""
            echo -e "${GREEN}✓ Ambiente actualizado exitosamente${NC}"
            echo ""
            echo "Herramientas instaladas:"
            micromamba list | grep -E 'metabat|maxbin|concoct|das_tool|bowtie2|samtools|bbmap|diamond'
            echo ""
            exit 0
        else
            echo ""
            echo -e "${RED}✗ Error al actualizar el ambiente${NC}"
            exit 1
        fi
    fi
fi

# Crear nuevo ambiente
echo -e "${CYAN}Creando nuevo ambiente 'binning'...${NC}"
echo ""

# Verificar si existe el archivo YAML
YAML_FILE="04_binning.yaml"

if [ -f "${YAML_FILE}" ]; then
    echo -e "${CYAN}Usando archivo YAML: ${YAML_FILE}${NC}"
    micromamba env create -f "${YAML_FILE}" -y
else
    echo -e "${YELLOW}⚠ No se encontró ${YAML_FILE}, creando ambiente manualmente${NC}"
    
    # Crear ambiente desde cero
    micromamba create -n binning -y -c bioconda -c conda-forge \
        bowtie2 \
        samtools \
        metabat2 \
        maxbin2 \
        concoct \
        das_tool \
        bbmap \
        diamond \
        checkm-genome
fi

if [ $? -eq 0 ]; then
    echo ""
    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}✓ Ambiente 'binning' creado exitosamente${NC}"
    echo -e "${GREEN}========================================${NC}"
    echo ""
    
    # Activar ambiente para verificar
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate binning
    
    echo "Herramientas instaladas:"
    echo ""
    
    # Verificar cada herramienta
    for tool in bowtie2-build bowtie2 samtools jgi_summarize_bam_contig_depths metabat2 run_MaxBin.pl cut_up_fasta.py concoct DAS_Tool pileup.sh diamond; do
        if command -v $tool &> /dev/null; then
            echo -e "  ${GREEN}✓${NC} $tool"
        else
            echo -e "  ${RED}✗${NC} $tool"
        fi
    done
    
    echo ""
    echo -e "${CYAN}Siguiente paso:${NC}"
    echo "  1. Activa el ambiente: ${YELLOW}micromamba activate binning${NC}"
    echo "  2. Ejecuta el binning: ${YELLOW}bash run_binning_fixed.sh${NC}"
    echo ""
else
    echo ""
    echo -e "${RED}✗ Error al crear el ambiente${NC}"
    echo ""
    echo "Solución:"
    echo "  1. Verifica tu conexión a internet"
    echo "  2. Intenta manualmente: micromamba create -n binning -c bioconda metabat2 maxbin2"
    echo ""
    exit 1
fi
