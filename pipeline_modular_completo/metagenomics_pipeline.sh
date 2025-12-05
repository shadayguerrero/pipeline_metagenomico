#!/bin/bash

# ============================================================================
# Pipeline Metagenómico Modular
# Autor: Shaday Guerrero / Manus AI
# Fecha: 2025-11-19
# Descripción: Pipeline completo para análisis metagenómico con selección
#              de módulos y ejecución secuencial
# ============================================================================

# ============================================================================
# COLORES PARA OUTPUT
# ============================================================================
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
NC='\033[0m'

# ============================================================================
# CONFIGURACIÓN GLOBAL
# ============================================================================

# Directorio base del proyecto
PROJECT_DIR="/files/shaday/4_cienegas"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Directorios de entrada/salida
INPUT_DIR="${PROJECT_DIR}/MergedFastq"
OUTPUT_DIR="${PROJECT_DIR}/output"

# Ambientes de micromamba (individuales)
ENV_BASE="/home_local/camda/micromamba/envs"
QC_ENV="${ENV_BASE}/qc_assembly"
BINNING_ENV="${ENV_BASE}/binning"
KRAKEN_ENV="kraken2"
PROKKA_ENV="prokka"
RGI_ENV="rgi"
ANTISMASH_ENV="antismash"
ANALYSIS_ENV="analysis"

# Bases de datos
BOWTIE2_INDEX="/home_local/camda/shaday/database/bowtie2/index/All/all"
KRAKEN2_GTDB="/files/database/k2_gtdb_r214"
KRAKEN2_PLUSPFP="/files/database/k2_pluspfp_20250402"
KRAKEN2_EUPATH="/files/database/k2_eupathdb48_20230407"

# Recursos computacionales
THREADS_QC=20
THREADS_HOST=12
THREADS_ASSEMBLY=60
THREADS_BINNING=40
THREADS_KRAKEN=40

# Directorios temporales
# IMPORTANTE: Si root está lleno, usar /mnt/Part4/Laboratory/tmp
# Ejecuta: source setup_tmp_part4.sh antes del pipeline
if [ -z "${TMPDIR}" ]; then
    # Intentar usar /mnt/Part4/Laboratory/tmp si existe
    if [ -d "/mnt/Part4/Laboratory" ]; then
        export TMPDIR="/mnt/Part4/Laboratory/tmp/general"
        mkdir -p "${TMPDIR}"
    elif [ -d "/Part4" ]; then
        export TMPDIR="/Part4/tmp/general"
        mkdir -p "${TMPDIR}"
    else
        export TMPDIR="/tmp"
    fi
fi
export TEMP="${TMPDIR}"
export TMP="${TMPDIR}"

# ============================================================================
# VARIABLES DE ESTADO
# ============================================================================
declare -A MODULES_SELECTED
declare -A MODULES_COMPLETED
KRAKEN_MODE=""  # "simple", "dual", "triple"
BINS_SOURCE="auto"  # "auto", "dastool", "metabat2", "maxbin2", "concoct"

# ============================================================================
# FUNCIONES DE UTILIDAD
# ============================================================================

print_header() {
    echo ""
    echo -e "${CYAN}============================================================================${NC}"
    echo -e "${CYAN}  $1${NC}"
    echo -e "${CYAN}============================================================================${NC}"
    echo ""
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}"
}

print_info() {
    echo -e "${BLUE}ℹ $1${NC}"
}

activate_env() {
    local env_name=$1
    print_info "Activando ambiente: ${env_name}"
    
    # Inicializar micromamba si no está ya inicializado
    if ! command -v micromamba &> /dev/null; then
        # Buscar micromamba en ubicaciones comunes
        if [ -f "${HOME}/micromamba/bin/micromamba" ]; then
            export PATH="${HOME}/micromamba/bin:$PATH"
        elif [ -f "/home_local/camda/micromamba/bin/micromamba" ]; then
            export PATH="/home_local/camda/micromamba/bin:$PATH"
        else
            print_error "No se encontró micromamba. Por favor, instálalo o añádelo al PATH."
            return 1
        fi
    fi
    
    # Inicializar el shell hook
    eval "$(micromamba shell hook --shell bash)"
    
    # Activar el ambiente
    micromamba activate "${env_name}" 2>/dev/null
    if [ $? -ne 0 ]; then
        print_warning "No se pudo activar el ambiente ${env_name}"
        print_info "Continuando sin ambiente activado..."
        return 1
    fi
    print_success "Ambiente ${env_name} activado"
}

deactivate_env() {
    if command -v micromamba &> /dev/null; then
        micromamba deactivate 2>/dev/null || true
    fi
}

check_file_exists() {
    if [ ! -f "$1" ]; then
        print_error "Archivo no encontrado: $1"
        return 1
    fi
    return 0
}

check_dir_exists() {
    if [ ! -d "$1" ]; then
        print_error "Directorio no encontrado: $1"
        return 1
    fi
    return 0
}

# ============================================================================
# MENÚ PRINCIPAL
# ============================================================================

show_main_menu() {
    clear
    print_header "PIPELINE METAGENÓMICO MODULAR"
    
    echo -e "${MAGENTA}Seleccione los módulos a ejecutar:${NC}"
    echo ""
    echo -e "  ${YELLOW}1.${NC} QC & Trimming (Trim Galore)"
    echo -e "  ${YELLOW}2.${NC} Host Removal (Bowtie2)"
    echo -e "  ${YELLOW}3.${NC} Assembly (MEGAHIT)"
    echo -e "  ${YELLOW}4.${NC} Binning (MetaBAT2/MaxBin2/CONCOCT)"
    echo -e "  ${YELLOW}5.${NC} Taxonomía de Reads (Kraken2)"
    echo -e "  ${YELLOW}6.${NC} Taxonomía de Bins (GTDB-Tk)"
    echo -e "  ${YELLOW}7.${NC} Anotación (Prokka)"
    echo -e "  ${YELLOW}8.${NC} Resistencia a Antibióticos (RGI)"
    echo -e "  ${YELLOW}9.${NC} Metabolitos Secundarios (AntiSMASH)"
    echo -e "  ${YELLOW}10.${NC} Análisis y Reportes"
    echo ""
    echo -e "  ${GREEN}A.${NC} Seleccionar TODOS los módulos"
    echo -e "  ${GREEN}B.${NC} Configurar fuente de bins (módulos 6-9)"
    echo -e "  ${GREEN}R.${NC} Revisar selección actual"
    echo -e "  ${GREEN}E.${NC} Ejecutar pipeline con módulos seleccionados"
    echo -e "  ${RED}Q.${NC} Salir"
    echo ""
}

show_kraken_menu() {
    clear
    print_header "CONFIGURACIÓN DE KRAKEN2"
    
    echo -e "${MAGENTA}Seleccione el modo de ejecución de Kraken2:${NC}"
    echo ""
    echo -e "  ${YELLOW}1.${NC} Simple - Solo GTDB (bacterias y arqueas)"
    echo -e "  ${YELLOW}2.${NC} Dual - GTDB + PlusPFP (+ eucariotas y virus)"
    echo -e "  ${YELLOW}3.${NC} Triple - GTDB + PlusPFP + EuPathDB (máxima cobertura)"
    echo ""
    echo -e "  ${RED}B.${NC} Volver al menú principal"
    echo ""
}

show_bins_menu() {
    clear
    print_header "SELECCIÓN DE BINS PARA ANÁLISIS"
    
    echo -e "${MAGENTA}Seleccione qué bins usar para módulos 6-9 (GTDB-Tk, Prokka, RGI, AntiSMASH):${NC}"
    echo ""
    echo -e "  ${YELLOW}1.${NC} Auto - Buscar automáticamente (DAS Tool → MetaBAT2 → MaxBin2 → CONCOCT)"
    echo -e "  ${YELLOW}2.${NC} DAS Tool - Solo bins refinados (RECOMENDADO para análisis final)"
    echo -e "  ${YELLOW}3.${NC} MetaBAT2 - Solo bins de MetaBAT2"
    echo -e "  ${YELLOW}4.${NC} MaxBin2 - Solo bins de MaxBin2"
    echo -e "  ${YELLOW}5.${NC} CONCOCT - Solo bins de CONCOCT"
    echo ""
    echo -e "${CYAN}Selección actual:${NC} ${BINS_SOURCE}"
    echo ""
    echo -e "  ${RED}B.${NC} Volver al menú principal"
    echo ""
}

toggle_module() {
    local module=$1
    if [ "${MODULES_SELECTED[$module]}" = "1" ]; then
        MODULES_SELECTED[$module]=0
        print_info "Módulo $module desactivado"
    else
        MODULES_SELECTED[$module]=1
        print_success "Módulo $module activado"
    fi
}

select_all_modules() {
    for i in {1..10}; do
        MODULES_SELECTED[$i]=1
    done
    print_success "Todos los módulos seleccionados"
}

show_selection() {
    print_header "MÓDULOS SELECCIONADOS"
    
    local modules=(
        "QC & Trimming"
        "Host Removal"
        "Assembly"
        "Binning"
        "Taxonomía Reads"
        "Taxonomía Bins"
        "Anotación"
        "RGI"
        "AntiSMASH"
        "Análisis"
    )
    
    for i in {1..10}; do
        if [ "${MODULES_SELECTED[$i]}" = "1" ]; then
            echo -e "  ${GREEN}✓${NC} ${modules[$((i-1))]}"
        else
            echo -e "  ${RED}✗${NC} ${modules[$((i-1))]}"
        fi
    done
    
    if [ "${MODULES_SELECTED[5]}" = "1" ]; then
        echo ""
        echo -e "${YELLOW}Modo Kraken2:${NC} ${KRAKEN_MODE:-No configurado}"
    fi
    
    # Mostrar fuente de bins si alguno de los módulos 6-9 está seleccionado
    if [ "${MODULES_SELECTED[6]}" = "1" ] || [ "${MODULES_SELECTED[7]}" = "1" ] || [ "${MODULES_SELECTED[8]}" = "1" ] || [ "${MODULES_SELECTED[9]}" = "1" ]; then
        echo ""
        echo -e "${CYAN}Fuente de bins (módulos 6-9):${NC} ${BINS_SOURCE}"
    fi
    
    echo ""
    read -p "Presione Enter para continuar..."
}

# ============================================================================
# MÓDULOS DEL PIPELINE
# ============================================================================

module_01_qc_trimming() {
    print_header "MÓDULO 1: QC & TRIMMING"
    
    activate_env "${QC_ENV}"
    
    local trim_output="${OUTPUT_DIR}/trim"
    mkdir -p "${trim_output}"
    
    print_info "Ejecutando Trim Galore..."
    
    cd "${INPUT_DIR}"
    for R1_FILE in *_1.fastq.gz; do
        SAMPLE=$(basename "${R1_FILE}" _1.fastq.gz)
        R2_FILE="${SAMPLE}_2.fastq.gz"
        
        if [ ! -f "${R2_FILE}" ]; then
            print_warning "No se encontró el archivo pareado ${R2_FILE}"
            continue
        fi
        
        print_info "Procesando muestra: ${SAMPLE}"
        
        trim_galore \
            --paired \
            --quality 20 \
            --length 20 \
            --cores ${THREADS_QC} \
            --fastqc \
            --output_dir "${trim_output}" \
            "${INPUT_DIR}/${R1_FILE}" \
            "${INPUT_DIR}/${R2_FILE}"
        
        if [ $? -eq 0 ]; then
            print_success "Muestra ${SAMPLE} completada"
        else
            print_error "Error en muestra ${SAMPLE}"
        fi
    done
    
    deactivate_env
    MODULES_COMPLETED[1]=1
    print_success "Módulo 1 completado"
}

module_02_host_removal() {
    print_header "MÓDULO 2: HOST REMOVAL"
    
    activate_env "${QC_ENV}"
    
    local input_dir="${OUTPUT_DIR}/trim"
    local output_dir="${OUTPUT_DIR}/host_removed"
    mkdir -p "${output_dir}"
    
    print_info "Removiendo reads del host con Bowtie2..."
    
    cd "${input_dir}"
    for R1_FILE in *_1_val_1.fq.gz; do
        SAMPLE=$(basename "${R1_FILE}" _1_val_1.fq.gz)
        R2_FILE="${SAMPLE}_2_val_2.fq.gz"
        
        if [ ! -f "${R2_FILE}" ]; then
            print_warning "No se encontró el archivo pareado ${R2_FILE}"
            continue
        fi
        
        print_info "Procesando muestra: ${SAMPLE}"
        
        local sample_output="${output_dir}/${SAMPLE}"
        mkdir -p "${sample_output}"
        
        # Mapeo
        bowtie2 -p ${THREADS_HOST} -x "${BOWTIE2_INDEX}" \
            -1 "${input_dir}/${R1_FILE}" \
            -2 "${input_dir}/${R2_FILE}" \
            -S "${sample_output}/${SAMPLE}_results.sam" \
            2> "${sample_output}/bowtie2.log"
        
        # SAM → BAM
        samtools view -bS "${sample_output}/${SAMPLE}_results.sam" \
            --threads ${THREADS_HOST} > "${sample_output}/${SAMPLE}_results.bam"
        
        # Filtrar no mapeados
        samtools view -b -f 13 -F 256 \
            "${sample_output}/${SAMPLE}_results.bam" \
            --threads ${THREADS_HOST} > "${sample_output}/${SAMPLE}_bothEndsUnmapped.bam"
        
        # Ordenar
        samtools sort -n -@ ${THREADS_HOST} \
            "${sample_output}/${SAMPLE}_bothEndsUnmapped.bam" \
            -o "${sample_output}/${SAMPLE}_sorted.bam"
        
        # BAM → FASTQ
        samtools fastq -@ ${THREADS_HOST} \
            "${sample_output}/${SAMPLE}_sorted.bam" \
            -1 "${sample_output}/${SAMPLE}_host_removed_R1.fastq.gz" \
            -2 "${sample_output}/${SAMPLE}_host_removed_R2.fastq.gz"
        
        # Limpiar archivos intermedios
        rm -f "${sample_output}/${SAMPLE}_results.sam" \
              "${sample_output}/${SAMPLE}_results.bam" \
              "${sample_output}/${SAMPLE}_bothEndsUnmapped.bam" \
              "${sample_output}/${SAMPLE}_sorted.bam"
        
        if [ $? -eq 0 ]; then
            print_success "Muestra ${SAMPLE} completada"
        else
            print_error "Error en muestra ${SAMPLE}"
        fi
    done
    
    deactivate_env
    MODULES_COMPLETED[2]=1
    print_success "Módulo 2 completado"
}

module_03_assembly() {
    print_header "MÓDULO 3: ASSEMBLY"
    
    activate_env "${QC_ENV}"
    
    local input_dir="${OUTPUT_DIR}/host_removed"
    local output_dir="${OUTPUT_DIR}/megahit_assemblies"
    mkdir -p "${output_dir}"
    
    print_info "Ensamblando con MEGAHIT..."
    
    for sample_dir in ${input_dir}/*; do
        if [ ! -d "${sample_dir}" ]; then
            continue
        fi
        
        SAMPLE=$(basename "${sample_dir}")
        R1_FILE="${sample_dir}/${SAMPLE}_host_removed_R1.fastq.gz"
        R2_FILE="${sample_dir}/${SAMPLE}_host_removed_R2.fastq.gz"
        
        if [ ! -f "${R1_FILE}" ] || [ ! -f "${R2_FILE}" ]; then
            print_warning "Archivos no encontrados para ${SAMPLE}"
            continue
        fi
        
        print_info "Ensamblando muestra: ${SAMPLE}"
        
        # Configurar directorio temporal para MEGAHIT
        MEGAHIT_TMP="${MEGAHIT_TMP:-${TMPDIR}/megahit}"
        mkdir -p "${MEGAHIT_TMP}"
        
        megahit \
            -1 "${R1_FILE}" \
            -2 "${R2_FILE}" \
            -o "${output_dir}/${SAMPLE}" \
            --tmp-dir "${MEGAHIT_TMP}" \
            --min-contig-len 500 \
            --k-min 21 \
            --k-max 141 \
            --k-step 12 \
            --preset meta-sensitive \
            -t ${THREADS_ASSEMBLY}
        
        if [ $? -eq 0 ]; then
            print_success "Muestra ${SAMPLE} completada"
        else
            print_error "Error en muestra ${SAMPLE}"
        fi
    done
    
    deactivate_env
    MODULES_COMPLETED[3]=1
    print_success "Módulo 3 completado"
}

module_04_binning() {
    print_header "MÓDULO 4: BINNING"
    
    activate_env "${BINNING_ENV}"
    
    local assembly_dir="${OUTPUT_DIR}/megahit_assemblies"
    local host_removed_dir="${OUTPUT_DIR}/host_removed"
    local output_dir="${OUTPUT_DIR}/binning"
    
    # Verificar directorios
    if [ ! -d "${assembly_dir}" ]; then
        print_error "No se encontró el directorio de ensamblajes: ${assembly_dir}"
        return 1
    fi
    
    if [ ! -d "${host_removed_dir}" ]; then
        print_error "No se encontró el directorio de reads: ${host_removed_dir}"
        return 1
    fi
    
    mkdir -p "${output_dir}"
    
    print_info "Ejecutando binning con MetaBAT2..."
    
    # Configurar variables de entorno para el script de binning
    export ASSEMBLY_DIR="${assembly_dir}"
    export HOST_REMOVED_DIR="${host_removed_dir}"
    export OUTPUT_DIR="${output_dir}"
    export THREADS="${THREADS_BINNING}"
    
    # Ejecutar script de binning
    if [ -f "${SCRIPT_DIR}/run_binning_fixed.sh" ]; then
        bash "${SCRIPT_DIR}/run_binning_fixed.sh"
    elif [ -f "${SCRIPT_DIR}/run_binning.sh" ]; then
        # Usar el script original si no existe el corregido
        bash "${SCRIPT_DIR}/run_binning.sh"
    else
        print_error "No se encontró el script de binning"
        return 1
    fi
    
    if [ $? -eq 0 ]; then
        print_success "Binning completado"
        MODULES_COMPLETED[4]=1
    else
        print_error "Error en binning"
        return 1
    fi
    
    deactivate_env
}

module_05_taxonomy_reads() {
    print_header "MÓDULO 5: TAXONOMÍA DE READS (KRAKEN2)"
    
    if [ -z "${KRAKEN_MODE}" ]; then
        print_error "Modo de Kraken2 no configurado"
        return 1
    fi
    
    activate_env "${KRAKEN_ENV}"
    
    local input_dir="${OUTPUT_DIR}/host_removed"
    local output_dir="${OUTPUT_DIR}/kraken2_${KRAKEN_MODE}"
    mkdir -p "${output_dir}"
    
    print_info "Ejecutando Kraken2 en modo: ${KRAKEN_MODE}"
    
    case "${KRAKEN_MODE}" in
        "simple")
            # Solo GTDB
            for sample_dir in ${input_dir}/*; do
                if [ ! -d "${sample_dir}" ]; then
                    continue
                fi
                
                SAMPLE=$(basename "${sample_dir}")
                R1_FILE="${sample_dir}/${SAMPLE}_host_removed_R1.fastq.gz"
                R2_FILE="${sample_dir}/${SAMPLE}_host_removed_R2.fastq.gz"
                
                print_info "Procesando muestra: ${SAMPLE}"
                
                kraken2 \
                    --db "${KRAKEN2_GTDB}" \
                    --paired "${R1_FILE}" "${R2_FILE}" \
                    --threads ${THREADS_KRAKEN} \
                    --report "${output_dir}/${SAMPLE}_GTDB.report" \
                    --output "${output_dir}/${SAMPLE}_GTDB.output"
                
                print_success "Muestra ${SAMPLE} completada"
            done
            ;;
        
        "dual")
            # GTDB + PlusPFP
            for sample_dir in ${input_dir}/*; do
                if [ ! -d "${sample_dir}" ]; then
                    continue
                fi
                
                SAMPLE=$(basename "${sample_dir}")
                R1_FILE="${sample_dir}/${SAMPLE}_host_removed_R1.fastq.gz"
                R2_FILE="${sample_dir}/${SAMPLE}_host_removed_R2.fastq.gz"
                
                print_info "Procesando muestra: ${SAMPLE} (GTDB)"
                
                # Paso 1: GTDB
                kraken2 \
                    --db "${KRAKEN2_GTDB}" \
                    --paired "${R1_FILE}" "${R2_FILE}" \
                    --threads ${THREADS_KRAKEN} \
                    --report "${output_dir}/${SAMPLE}_GTDB.report" \
                    --output "${output_dir}/${SAMPLE}_GTDB.output" \
                    --unclassified-out "${output_dir}/${SAMPLE}_GTDB_unclassified#.fastq"
                
                # Paso 2: PlusPFP con no clasificados
                print_info "Procesando muestra: ${SAMPLE} (PlusPFP)"
                
                kraken2 \
                    --db "${KRAKEN2_PLUSPFP}" \
                    --paired "${output_dir}/${SAMPLE}_GTDB_unclassified_1.fastq" \
                             "${output_dir}/${SAMPLE}_GTDB_unclassified_2.fastq" \
                    --threads ${THREADS_KRAKEN} \
                    --report "${output_dir}/${SAMPLE}_PlusPFP.report" \
                    --output "${output_dir}/${SAMPLE}_PlusPFP.output"
                
                # Limpiar archivos intermedios
                rm -f "${output_dir}/${SAMPLE}_GTDB_unclassified_"*.fastq
                
                print_success "Muestra ${SAMPLE} completada"
            done
            ;;
        
        "triple")
            # GTDB + PlusPFP + EuPathDB
            for sample_dir in ${input_dir}/*; do
                if [ ! -d "${sample_dir}" ]; then
                    continue
                fi
                
                SAMPLE=$(basename "${sample_dir}")
                R1_FILE="${sample_dir}/${SAMPLE}_host_removed_R1.fastq.gz"
                R2_FILE="${sample_dir}/${SAMPLE}_host_removed_R2.fastq.gz"
                
                print_info "Procesando muestra: ${SAMPLE} (GTDB)"
                
                # Paso 1: GTDB
                kraken2 \
                    --db "${KRAKEN2_GTDB}" \
                    --paired "${R1_FILE}" "${R2_FILE}" \
                    --threads ${THREADS_KRAKEN} \
                    --report "${output_dir}/${SAMPLE}_GTDB.report" \
                    --output "${output_dir}/${SAMPLE}_GTDB.output" \
                    --unclassified-out "${output_dir}/${SAMPLE}_GTDB_unclassified#.fastq"
                
                # Paso 2: PlusPFP
                print_info "Procesando muestra: ${SAMPLE} (PlusPFP)"
                
                kraken2 \
                    --db "${KRAKEN2_PLUSPFP}" \
                    --paired "${output_dir}/${SAMPLE}_GTDB_unclassified_1.fastq" \
                             "${output_dir}/${SAMPLE}_GTDB_unclassified_2.fastq" \
                    --threads ${THREADS_KRAKEN} \
                    --report "${output_dir}/${SAMPLE}_PlusPFP.report" \
                    --output "${output_dir}/${SAMPLE}_PlusPFP.output" \
                    --unclassified-out "${output_dir}/${SAMPLE}_PlusPFP_unclassified#.fastq"
                
                # Paso 3: EuPathDB
                print_info "Procesando muestra: ${SAMPLE} (EuPathDB)"
                
                kraken2 \
                    --db "${KRAKEN2_EUPATH}" \
                    --paired "${output_dir}/${SAMPLE}_PlusPFP_unclassified_1.fastq" \
                             "${output_dir}/${SAMPLE}_PlusPFP_unclassified_2.fastq" \
                    --threads ${THREADS_KRAKEN} \
                    --report "${output_dir}/${SAMPLE}_EuPathDB.report" \
                    --output "${output_dir}/${SAMPLE}_EuPathDB.output"
                
                # Limpiar archivos intermedios
                rm -f "${output_dir}/${SAMPLE}_GTDB_unclassified_"*.fastq
                rm -f "${output_dir}/${SAMPLE}_PlusPFP_unclassified_"*.fastq
                
                print_success "Muestra ${SAMPLE} completada"
            done
            ;;
    esac
    
    deactivate_env
    MODULES_COMPLETED[5]=1
    print_success "Módulo 5 completado"
}

module_06_gtdbtk() {
    print_header "MÓDULO 6: TAXONOMÍA DE BINS (GTDB-TK)"
    
    # Exportar fuente de bins
    export BINS_SOURCE
    
    # Configurar directorio temporal para GTDB-Tk
    export GTDBTK_TMP="${GTDBTK_TMP:-${TMPDIR}/gtdbtk}"
    mkdir -p "${GTDBTK_TMP}"
    
    activate_env "gtdbtk"
    
    local binning_dir="${OUTPUT_DIR}/binning"
    local output_dir="${OUTPUT_DIR}/gtdbtk"
    
    # Verificar directorio de binning
    if [ ! -d "${binning_dir}" ]; then
        print_error "No se encontró el directorio de binning: ${binning_dir}"
        print_info "Ejecuta primero el Módulo 4 (Binning)"
        return 1
    fi
    
    mkdir -p "${output_dir}"
    
    print_info "Clasificando bins con GTDB-Tk..."
    
    # Verificar que GTDB-Tk esté disponible
    if ! command -v gtdbtk &> /dev/null; then
        print_error "GTDB-Tk no está disponible"
        print_info "Instala: micromamba install -n gtdbtk gtdbtk"
        return 1
    fi
    
    # Verificar base de datos
    if [ -z "${GTDBTK_DATA_PATH}" ]; then
        print_error "Variable GTDBTK_DATA_PATH no configurada"
        print_info "Descarga la base de datos: https://data.gtdb.ecogenomic.org/"
        return 1
    fi
    
    # Procesar cada muestra
    for sample_dir in ${binning_dir}/*; do
        if [ ! -d "${sample_dir}" ]; then
            continue
        fi
        
        SAMPLE=$(basename "${sample_dir}")
        bins_dir="${sample_dir}/metabat2"
        
        if [ ! -d "${bins_dir}" ]; then
            print_warning "No se encontraron bins para ${SAMPLE}"
            continue
        fi
        
        # Contar bins
        bin_count=$(ls ${bins_dir}/bin.*.fa 2>/dev/null | wc -l)
        if [ ${bin_count} -eq 0 ]; then
            print_warning "No hay bins para ${SAMPLE}"
            continue
        fi
        
        print_info "Procesando ${SAMPLE} (${bin_count} bins)..."
        
        sample_output="${output_dir}/${SAMPLE}"
        mkdir -p "${sample_output}"
        
        # Ejecutar GTDB-Tk
        gtdbtk classify_wf \
            --genome_dir "${bins_dir}" \
            --out_dir "${sample_output}" \
            --extension fa \
            --cpus ${THREADS_BINNING} \
            --skip_ani_screen \
            > "${sample_output}/gtdbtk.log" 2>&1
        
        if [ $? -eq 0 ]; then
            print_success "Muestra ${SAMPLE} completada"
        else
            print_error "Error en muestra ${SAMPLE}"
        fi
    done
    
    deactivate_env
    MODULES_COMPLETED[6]=1
    print_success "Módulo 6 completado"
}

module_07_prokka() {
    print_header "MÓDULO 7: ANOTACIÓN (PROKKA)"
    
    # Exportar fuente de bins
    export BINS_SOURCE
    
    activate_env "prokka"
    
    local binning_dir="${OUTPUT_DIR}/binning"
    local gtdbtk_dir="${OUTPUT_DIR}/gtdbtk"
    local output_dir="${OUTPUT_DIR}/prokka"
    
    # Verificar directorios
    if [ ! -d "${binning_dir}" ]; then
        print_error "No se encontró el directorio de binning: ${binning_dir}"
        return 1
    fi
    
    if [ ! -d "${gtdbtk_dir}" ]; then
        print_warning "No se encontró GTDB-Tk, se usará reino por defecto (Bacteria)"
    fi
    
    mkdir -p "${output_dir}"
    
    print_info "Anotando bins con Prokka..."
    
    # Verificar que Prokka esté disponible
    if ! command -v prokka &> /dev/null; then
        print_error "Prokka no está disponible"
        return 1
    fi
    
    # Procesar cada muestra
    for sample_dir in ${binning_dir}/*; do
        if [ ! -d "${sample_dir}" ]; then
            continue
        fi
        
        SAMPLE=$(basename "${sample_dir}")
        bins_dir="${sample_dir}/metabat2"
        
        if [ ! -d "${bins_dir}" ]; then
            continue
        fi
        
        print_info "Procesando ${SAMPLE}..."
        
        sample_output="${output_dir}/${SAMPLE}"
        mkdir -p "${sample_output}"
        
        # Procesar cada bin
        for bin_file in ${bins_dir}/bin.*.fa; do
            if [ ! -f "${bin_file}" ]; then
                continue
            fi
            
            bin_name=$(basename "${bin_file}" .fa)
            bin_output="${sample_output}/${bin_name}"
            mkdir -p "${bin_output}"
            
            # Determinar reino (Bacteria por defecto)
            kingdom="Bacteria"
            
            # Buscar clasificación en GTDB-Tk si existe
            if [ -f "${gtdbtk_dir}/${SAMPLE}/gtdbtk.bac120.summary.tsv" ]; then
                kingdom="Bacteria"
            elif [ -f "${gtdbtk_dir}/${SAMPLE}/gtdbtk.ar53.summary.tsv" ] || [ -f "${gtdbtk_dir}/${SAMPLE}/gtdbtk.ar122.summary.tsv" ]; then
                kingdom="Archaea"
            fi
            
            print_info "  Anotando ${bin_name} (${kingdom})..."
            
            prokka \
                --outdir "${bin_output}" \
                --prefix "${bin_name}" \
                --cpus 4 \
                --mincontiglen 200 \
                --metagenome \
                --force \
                --kingdom "${kingdom}" \
                "${bin_file}" > "${bin_output}/prokka.log" 2>&1
            
            if [ $? -eq 0 ]; then
                print_success "  ${bin_name} completado"
            else
                print_warning "  Error en ${bin_name}"
            fi
        done
    done
    
    deactivate_env
    MODULES_COMPLETED[7]=1
    print_success "Módulo 7 completado"
}

module_08_rgi() {
    print_header "MÓDULO 8: RESISTENCIA A ANTIBIÓTICOS (RGI)"
    
    # Exportar fuente de bins
    export BINS_SOURCE
    
    activate_env "rgi"
    
    local prokka_dir="${OUTPUT_DIR}/prokka"
    local output_dir="${OUTPUT_DIR}/rgi"
    
    # Verificar directorio de Prokka
    if [ ! -d "${prokka_dir}" ]; then
        print_error "No se encontró el directorio de Prokka: ${prokka_dir}"
        print_info "Ejecuta primero el Módulo 7 (Prokka)"
        return 1
    fi
    
    mkdir -p "${output_dir}"
    
    print_info "Analizando resistencia con RGI..."
    
    # Verificar que RGI esté disponible
    if ! command -v rgi &> /dev/null; then
        print_error "RGI no está disponible"
        return 1
    fi
    
    # Verificar base de datos CARD
    print_info "Verificando base de datos CARD..."
    rgi load --help > /dev/null 2>&1 || {
        print_warning "Base de datos CARD no cargada, intentando cargar..."
        rgi load --card_json /path/to/card.json || {
            print_error "No se pudo cargar la base de datos CARD"
            print_info "Descarga desde: https://card.mcmaster.ca/download"
            return 1
        }
    }
    
    # Procesar cada muestra
    for sample_dir in ${prokka_dir}/*; do
        if [ ! -d "${sample_dir}" ]; then
            continue
        fi
        
        SAMPLE=$(basename "${sample_dir}")
        print_info "Procesando ${SAMPLE}..."
        
        sample_output="${output_dir}/${SAMPLE}"
        mkdir -p "${sample_output}"
        
        # Procesar cada bin
        for bin_dir in ${sample_dir}/bin.*; do
            if [ ! -d "${bin_dir}" ]; then
                continue
            fi
            
            bin_name=$(basename "${bin_dir}")
            faa_file="${bin_dir}/${bin_name}.faa"
            
            if [ ! -f "${faa_file}" ]; then
                continue
            fi
            
            print_info "  Analizando ${bin_name}..."
            
            bin_output="${sample_output}/${bin_name}"
            mkdir -p "${bin_output}"
            
            rgi main \
                --input_sequence "${faa_file}" \
                --output_file "${bin_output}/${bin_name}" \
                --input_type protein \
                --num_threads 4 \
                --clean \
                > "${bin_output}/rgi.log" 2>&1
            
            if [ $? -eq 0 ]; then
                print_success "  ${bin_name} completado"
            else
                print_warning "  Error en ${bin_name}"
            fi
        done
    done
    
    deactivate_env
    MODULES_COMPLETED[8]=1
    print_success "Módulo 8 completado"
}

module_09_antismash() {
    print_header "MÓDULO 9: METABOLITOS SECUNDARIOS (ANTISMASH)"
    
    # Exportar fuente de bins
    export BINS_SOURCE
    
    activate_env "antismash"
    
    local prokka_dir="${OUTPUT_DIR}/prokka"
    local output_dir="${OUTPUT_DIR}/antismash"
    
    # Verificar directorio de Prokka
    if [ ! -d "${prokka_dir}" ]; then
        print_error "No se encontró el directorio de Prokka: ${prokka_dir}"
        print_info "Ejecuta primero el Módulo 7 (Prokka)"
        return 1
    fi
    
    mkdir -p "${output_dir}"
    
    print_info "Analizando metabolitos secundarios con AntiSMASH..."
    
    # Verificar que AntiSMASH esté disponible
    if ! command -v antismash &> /dev/null; then
        print_error "AntiSMASH no está disponible"
        return 1
    fi
    
    # Procesar cada muestra
    for sample_dir in ${prokka_dir}/*; do
        if [ ! -d "${sample_dir}" ]; then
            continue
        fi
        
        SAMPLE=$(basename "${sample_dir}")
        print_info "Procesando ${SAMPLE}..."
        
        sample_output="${output_dir}/${SAMPLE}"
        mkdir -p "${sample_output}"
        
        # Procesar cada bin
        for bin_dir in ${sample_dir}/bin.*; do
            if [ ! -d "${bin_dir}" ]; then
                continue
            fi
            
            bin_name=$(basename "${bin_dir}")
            gbk_file="${bin_dir}/${bin_name}.gbk"
            
            if [ ! -f "${gbk_file}" ]; then
                continue
            fi
            
            print_info "  Analizando ${bin_name}..."
            
            bin_output="${sample_output}/${bin_name}"
            
            antismash \
                --genefinding-tool prodigal \
                --output-dir "${bin_output}" \
                --cpus 4 \
                --taxon bacteria \
                --minlength 1000 \
                "${gbk_file}" > "${bin_output}/antismash.log" 2>&1
            
            if [ $? -eq 0 ]; then
                print_success "  ${bin_name} completado"
            else
                print_warning "  Error en ${bin_name}"
            fi
        done
    done
    
    deactivate_env
    MODULES_COMPLETED[9]=1
    print_success "Módulo 9 completado"
}

module_10_analysis() {
    print_header "MÓDULO 10: ANÁLISIS Y REPORTES"
    
    activate_env "${ANALYSIS_ENV}"
    
    local kraken_dir="${OUTPUT_DIR}/kraken2_${KRAKEN_MODE}"
    local output_dir="${OUTPUT_DIR}/analysis"
    mkdir -p "${output_dir}"
    
    print_info "Generando archivo BIOM..."
    
    case "${KRAKEN_MODE}" in
        "simple")
            python3 "${SCRIPT_DIR}/combinar_kraken_simple_a_biom.py" \
                "${kraken_dir}" \
                "${output_dir}/combined_results.biom"
            ;;
        "dual")
            python3 "${SCRIPT_DIR}/combinar_kraken_2bases_a_biom.py" \
                "${kraken_dir}" \
                "${output_dir}/combined_results.biom"
            ;;
        "triple")
            python3 "${SCRIPT_DIR}/combinar_kraken_3bases_a_biom_v4.py" \
                "${kraken_dir}" \
                "${output_dir}/combined_results.biom"
            ;;
    esac
    
    if [ $? -ne 0 ]; then
        print_error "Error al generar archivo BIOM"
        deactivate_env
        return 1
    fi
    
    print_success "Archivo BIOM generado"
    
    print_info "Generando reporte HTML..."
    
    python3 "${SCRIPT_DIR}/analisis_metagenomico_completo.py" \
        "${output_dir}/combined_results.biom" \
        "${PROJECT_DIR}/metadata.csv"
    
    if [ $? -eq 0 ]; then
        print_success "Reporte HTML generado"
        print_info "Reporte disponible en: ${output_dir}/resultados_metagenomicos/reporte_metagenomico.html"
    else
        print_error "Error al generar reporte HTML"
    fi
    
    deactivate_env
    MODULES_COMPLETED[10]=1
    print_success "Módulo 10 completado"
}

# ============================================================================
# EJECUCIÓN DEL PIPELINE
# ============================================================================

execute_pipeline() {
    print_header "EJECUTANDO PIPELINE"
    
    # Verificar que al menos un módulo esté seleccionado
    local has_selection=0
    for i in {1..10}; do
        if [ "${MODULES_SELECTED[$i]}" = "1" ]; then
            has_selection=1
            break
        fi
    done
    
    if [ $has_selection -eq 0 ]; then
        print_error "No hay módulos seleccionados"
        return 1
    fi
    
    # Verificar configuración de Kraken2 si está seleccionado
    if [ "${MODULES_SELECTED[5]}" = "1" ] && [ -z "${KRAKEN_MODE}" ]; then
        print_error "Kraken2 seleccionado pero no configurado"
        return 1
    fi
    
    # Crear directorio de salida
    mkdir -p "${OUTPUT_DIR}"
    
    # Ejecutar módulos en orden
    local start_time=$(date +%s)
    
    [ "${MODULES_SELECTED[1]}" = "1" ] && module_01_qc_trimming
    [ "${MODULES_SELECTED[2]}" = "1" ] && module_02_host_removal
    [ "${MODULES_SELECTED[3]}" = "1" ] && module_03_assembly
    [ "${MODULES_SELECTED[4]}" = "1" ] && module_04_binning
    [ "${MODULES_SELECTED[5]}" = "1" ] && module_05_taxonomy_reads
    [ "${MODULES_SELECTED[6]}" = "1" ] && module_06_gtdbtk
    [ "${MODULES_SELECTED[7]}" = "1" ] && module_07_prokka
    [ "${MODULES_SELECTED[8]}" = "1" ] && module_08_rgi
    [ "${MODULES_SELECTED[9]}" = "1" ] && module_09_antismash
    [ "${MODULES_SELECTED[10]}" = "1" ] && module_10_analysis
    
    local end_time=$(date +%s)
    local elapsed=$((end_time - start_time))
    local hours=$((elapsed / 3600))
    local minutes=$(((elapsed % 3600) / 60))
    local seconds=$((elapsed % 60))
    
    print_header "PIPELINE COMPLETADO"
    print_success "Tiempo total: ${hours}h ${minutes}m ${seconds}s"
    
    echo ""
    read -p "Presione Enter para volver al menú principal..."
}

# ============================================================================
# BUCLE PRINCIPAL
# ============================================================================

main() {
    # Inicializar variables
    for i in {1..10}; do
        MODULES_SELECTED[$i]=0
        MODULES_COMPLETED[$i]=0
    done
    
    while true; do
        show_main_menu
        read -p "Seleccione una opción: " option
        
        case $option in
            [1-9])
                toggle_module $option
                ;;
            10)
                toggle_module 10
                ;;
            A|a)
                select_all_modules
                sleep 1
                ;;
            B|b)
                # Configurar fuente de bins
                while true; do
                    show_bins_menu
                    read -p "Seleccione una opción: " bins_option
                    
                    case $bins_option in
                        1)
                            BINS_SOURCE="auto"
                            print_success "Fuente de bins: Auto (búsqueda automática)"
                            sleep 1
                            ;;
                        2)
                            BINS_SOURCE="dastool"
                            print_success "Fuente de bins: DAS Tool (refinados)"
                            sleep 1
                            ;;
                        3)
                            BINS_SOURCE="metabat2"
                            print_success "Fuente de bins: MetaBAT2"
                            sleep 1
                            ;;
                        4)
                            BINS_SOURCE="maxbin2"
                            print_success "Fuente de bins: MaxBin2"
                            sleep 1
                            ;;
                        5)
                            BINS_SOURCE="concoct"
                            print_success "Fuente de bins: CONCOCT"
                            sleep 1
                            ;;
                        B|b)
                            break
                            ;;
                        *)
                            print_error "Opción inválida"
                            sleep 1
                            ;;
                    esac
                done
                ;;
            R|r)
                show_selection
                ;;
            E|e)
                # Si Kraken2 está seleccionado y no configurado, mostrar menú
                if [ "${MODULES_SELECTED[5]}" = "1" ] && [ -z "${KRAKEN_MODE}" ]; then
                    while true; do
                        show_kraken_menu
                        read -p "Seleccione una opción: " kraken_option
                        
                        case $kraken_option in
                            1)
                                KRAKEN_MODE="simple"
                                print_success "Modo Kraken2: Simple (GTDB)"
                                sleep 1
                                break
                                ;;
                            2)
                                KRAKEN_MODE="dual"
                                print_success "Modo Kraken2: Dual (GTDB + PlusPFP)"
                                sleep 1
                                break
                                ;;
                            3)
                                KRAKEN_MODE="triple"
                                print_success "Modo Kraken2: Triple (GTDB + PlusPFP + EuPathDB)"
                                sleep 1
                                break
                                ;;
                            B|b)
                                break
                                ;;
                            *)
                                print_error "Opción inválida"
                                sleep 1
                                ;;
                        esac
                    done
                fi
                
                execute_pipeline
                ;;
            Q|q)
                print_info "Saliendo del pipeline..."
                exit 0
                ;;
            *)
                print_error "Opción inválida"
                sleep 1
                ;;
        esac
    done
}

# ============================================================================
# PUNTO DE ENTRADA
# ============================================================================

main "$@"

