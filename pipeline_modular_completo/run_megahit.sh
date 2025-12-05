#!/bin/bash

# Script para ensamblaje metagenómico con MEGAHIT
# Procesa archivos FASTQ pareados sin contaminación del host
# Autor: Script generado para análisis de 4 Ciénegas
# Fecha: $(date +%Y-%m-%d)

# ============================================================================
# CONFIGURACIÓN
# ============================================================================

# Directorios
INPUT_DIR="/files/shaday/4_cienegas/output/host_removed/"
OUTPUT_DIR="/files/shaday/4_cienegas/output/megahit_assemblies"

# Parámetros de procesamiento
THREADS_PER_SAMPLE=60  # Threads por muestra (MEGAHIT usa mucha memoria y CPU)
MAX_PARALLEL_JOBS=1    # Número de muestras en paralelo (recomendado: 1-2 por limitaciones de memoria)

# Parámetros de MEGAHIT
MIN_CONTIG_LEN=500     # Longitud mínima de contigs (bp)
K_MIN=21               # K-mer mínimo
K_MAX=141              # K-mer máximo
K_STEP=12              # Incremento de k-mer

# Preset de MEGAHIT (meta-sensitive, meta-large, o vacío para default)
PRESET="meta-sensitive"  # Opciones: "", "meta-sensitive", "meta-large"

# Sufijos de archivos de host removal
SUFFIX_R1="_host_removed_R1.fastq.gz"
SUFFIX_R2="_host_removed_R2.fastq.gz"

# Colores para output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
NC='\033[0m' # No Color

# ============================================================================
# FUNCIONES
# ============================================================================

# Función para procesar una muestra individual con MEGAHIT
process_megahit() {
    local R1_FILE=$1
    local SAMPLE=$2
    local R2_FILE=$3
    local THREADS=$4
    local SAMPLE_OUTPUT="${OUTPUT_DIR}/${SAMPLE}"
    
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}[$(date +%H:%M:%S)] Procesando: ${YELLOW}${SAMPLE}${NC}"
    echo -e "${BLUE}========================================${NC}"
    
    # Verificar si ya existe un ensamblaje previo
    if [ -d "${SAMPLE_OUTPUT}" ]; then
        echo -e "${YELLOW}  ⚠ Directorio de salida ya existe: ${SAMPLE_OUTPUT}${NC}"
        echo -e "${YELLOW}  Se eliminará y se reiniciará el ensamblaje${NC}"
        rm -rf "${SAMPLE_OUTPUT}"
    fi
    
    # ========================================================================
    # ENSAMBLAJE CON MEGAHIT
    # ========================================================================
    echo -e "${CYAN}  [1/3] Iniciando ensamblaje con MEGAHIT...${NC}"
    echo -e "${CYAN}    R1: $(basename ${R1_FILE})${NC}"
    echo -e "${CYAN}    R2: $(basename ${R2_FILE})${NC}"
    echo -e "${CYAN}    Threads: ${THREADS}${NC}"
    echo -e "${CYAN}    K-mer range: ${K_MIN}-${K_MAX} (step ${K_STEP})${NC}"
    echo -e "${CYAN}    Min contig length: ${MIN_CONTIG_LEN} bp${NC}"
    if [ -n "${PRESET}" ]; then
        echo -e "${CYAN}    Preset: ${PRESET}${NC}"
    fi
    echo ""
    
    START_TIME=$(date +%s)
    
    # Construir comando MEGAHIT
    MEGAHIT_CMD="megahit \
        -1 ${R1_FILE} \
        -2 ${R2_FILE} \
        -o ${SAMPLE_OUTPUT} \
        -t ${THREADS} \
        --min-contig-len ${MIN_CONTIG_LEN} \
        --k-min ${K_MIN} \
        --k-max ${K_MAX} \
        --k-step ${K_STEP}"
    
    # Agregar preset si está definido
    if [ -n "${PRESET}" ]; then
        MEGAHIT_CMD="${MEGAHIT_CMD} --presets ${PRESET}"
    fi
    
    # Ejecutar MEGAHIT
    eval ${MEGAHIT_CMD} 2>&1 | tee "${OUTPUT_DIR}/${SAMPLE}_megahit.log"
    
    MEGAHIT_EXIT=${PIPESTATUS[0]}
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    ELAPSED_HOURS=$((ELAPSED / 3600))
    ELAPSED_MIN=$(((ELAPSED % 3600) / 60))
    ELAPSED_SEC=$((ELAPSED % 60))
    
    if [ ${MEGAHIT_EXIT} -ne 0 ]; then
        echo -e "${RED}  ✗ Error en MEGAHIT para ${SAMPLE}${NC}"
        echo -e "${RED}    Ver log: ${OUTPUT_DIR}/${SAMPLE}_megahit.log${NC}"
        return 1
    fi
    
    echo ""
    echo -e "${GREEN}  ✓ Ensamblaje completado (${ELAPSED_HOURS}h ${ELAPSED_MIN}m ${ELAPSED_SEC}s)${NC}"
    
    # ========================================================================
    # ANÁLISIS DE ESTADÍSTICAS DEL ENSAMBLAJE
    # ========================================================================
    echo -e "${CYAN}  [2/3] Calculando estadísticas del ensamblaje...${NC}"
    
    FINAL_CONTIGS="${SAMPLE_OUTPUT}/final.contigs.fa"
    
    if [ ! -f "${FINAL_CONTIGS}" ]; then
        echo -e "${RED}  ✗ No se encontró el archivo de contigs: ${FINAL_CONTIGS}${NC}"
        return 1
    fi
    
    # Calcular estadísticas con script Python
    python3 << PYTHON_STATS
import sys
from collections import defaultdict

def calculate_assembly_stats(fasta_file):
    lengths = []
    total_bases = 0
    gc_count = 0
    
    with open(fasta_file, 'r') as f:
        seq = ""
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq:
                    lengths.append(len(seq))
                    total_bases += len(seq)
                    gc_count += seq.upper().count('G') + seq.upper().count('C')
                seq = ""
            else:
                seq += line
        # Última secuencia
        if seq:
            lengths.append(len(seq))
            total_bases += len(seq)
            gc_count += seq.upper().count('G') + seq.upper().count('C')
    
    if not lengths:
        return None
    
    # Ordenar longitudes de mayor a menor
    lengths.sort(reverse=True)
    
    # Estadísticas básicas
    num_contigs = len(lengths)
    total_length = sum(lengths)
    longest = lengths[0]
    shortest = lengths[-1]
    mean_length = total_length / num_contigs
    gc_content = (gc_count / total_bases * 100) if total_bases > 0 else 0
    
    # Calcular N50, N75, N90
    def calculate_nx(lengths, x):
        target = sum(lengths) * x / 100
        cumsum = 0
        for length in lengths:
            cumsum += length
            if cumsum >= target:
                return length
        return lengths[-1]
    
    n50 = calculate_nx(lengths, 50)
    n75 = calculate_nx(lengths, 75)
    n90 = calculate_nx(lengths, 90)
    
    # Calcular L50 (número de contigs que suman al 50% del total)
    l50 = 0
    cumsum = 0
    target = total_length * 0.5
    for length in lengths:
        l50 += 1
        cumsum += length
        if cumsum >= target:
            break
    
    # Contigs por rango de tamaño
    ranges = {
        '500-1000': 0,
        '1000-5000': 0,
        '5000-10000': 0,
        '10000-50000': 0,
        '>50000': 0
    }
    
    for length in lengths:
        if 500 <= length < 1000:
            ranges['500-1000'] += 1
        elif 1000 <= length < 5000:
            ranges['1000-5000'] += 1
        elif 5000 <= length < 10000:
            ranges['5000-10000'] += 1
        elif 10000 <= length < 50000:
            ranges['10000-50000'] += 1
        elif length >= 50000:
            ranges['>50000'] += 1
    
    return {
        'num_contigs': num_contigs,
        'total_length': total_length,
        'longest': longest,
        'shortest': shortest,
        'mean_length': mean_length,
        'n50': n50,
        'n75': n75,
        'n90': n90,
        'l50': l50,
        'gc_content': gc_content,
        'ranges': ranges
    }

# Calcular estadísticas
stats = calculate_assembly_stats("${FINAL_CONTIGS}")

if stats:
    # Guardar estadísticas en archivo
    with open("${SAMPLE_OUTPUT}/assembly_stats.txt", 'w') as out:
        out.write("=" * 60 + "\n")
        out.write(f"Estadísticas del Ensamblaje - ${SAMPLE}\n")
        out.write("=" * 60 + "\n")
        out.write(f"Fecha: $(date +%Y-%m-%d\ %H:%M:%S)\n")
        out.write(f"Tiempo de ensamblaje: ${ELAPSED_HOURS}h ${ELAPSED_MIN}m ${ELAPSED_SEC}s\n\n")
        
        out.write("ESTADÍSTICAS GENERALES:\n")
        out.write(f"  Número de contigs: {stats['num_contigs']:,}\n")
        out.write(f"  Longitud total: {stats['total_length']:,} bp ({stats['total_length']/1e6:.2f} Mbp)\n")
        out.write(f"  Longitud promedio: {stats['mean_length']:.0f} bp\n")
        out.write(f"  Contig más largo: {stats['longest']:,} bp\n")
        out.write(f"  Contig más corto: {stats['shortest']:,} bp\n")
        out.write(f"  Contenido GC: {stats['gc_content']:.2f}%\n\n")
        
        out.write("MÉTRICAS DE CALIDAD:\n")
        out.write(f"  N50: {stats['n50']:,} bp\n")
        out.write(f"  N75: {stats['n75']:,} bp\n")
        out.write(f"  N90: {stats['n90']:,} bp\n")
        out.write(f"  L50: {stats['l50']:,} contigs\n\n")
        
        out.write("DISTRIBUCIÓN POR TAMAÑO:\n")
        for range_name, count in stats['ranges'].items():
            out.write(f"  {range_name} bp: {count:,} contigs\n")
        
        out.write("\n" + "=" * 60 + "\n")
        out.write("PARÁMETROS DE MEGAHIT:\n")
        out.write(f"  K-mer range: ${K_MIN}-${K_MAX} (step ${K_STEP})\n")
        out.write(f"  Min contig length: ${MIN_CONTIG_LEN} bp\n")
        out.write(f"  Preset: ${PRESET}\n")
        out.write(f"  Threads: ${THREADS}\n")
        out.write("\n" + "=" * 60 + "\n")
    
    # Imprimir resumen en consola
    print(f"\n{stats['num_contigs']:,} contigs, {stats['total_length']/1e6:.2f} Mbp total")
    print(f"N50: {stats['n50']:,} bp, L50: {stats['l50']:,} contigs")
    print(f"Longest: {stats['longest']:,} bp, GC: {stats['gc_content']:.2f}%")
else:
    print("Error al calcular estadísticas")
    sys.exit(1)
PYTHON_STATS
    
    STATS_EXIT=$?
    
    if [ ${STATS_EXIT} -ne 0 ]; then
        echo -e "${YELLOW}  ⚠ No se pudieron calcular estadísticas completas${NC}"
    else
        echo -e "${GREEN}  ✓ Estadísticas calculadas${NC}"
    fi
    
    # ========================================================================
    # LIMPIEZA DE ARCHIVOS INTERMEDIOS (OPCIONAL)
    # ========================================================================
    echo -e "${CYAN}  [3/3] Limpieza de archivos intermedios...${NC}"
    
    # MEGAHIT genera muchos archivos intermedios que pueden eliminarse
    # Descomentar las siguientes líneas si necesitas ahorrar espacio
    
    # rm -rf "${SAMPLE_OUTPUT}/intermediate_contigs"
    # rm -f "${SAMPLE_OUTPUT}/*.tmp"
    
    echo -e "${GREEN}  ✓ Limpieza completada${NC}"
    
    # ========================================================================
    # RESUMEN FINAL
    # ========================================================================
    echo ""
    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}  ✓ ${SAMPLE} completado exitosamente${NC}"
    echo -e "${GREEN}========================================${NC}"
    
    if [ -f "${SAMPLE_OUTPUT}/assembly_stats.txt" ]; then
        echo -e "${CYAN}Resumen:${NC}"
        grep -E "Número de contigs|Longitud total|N50|Contig más largo" "${SAMPLE_OUTPUT}/assembly_stats.txt" | sed 's/^/  /'
    fi
    
    echo ""
    
    return 0
}

# Exportar función y variables para GNU Parallel
export -f process_megahit
export OUTPUT_DIR THREADS_PER_SAMPLE
export MIN_CONTIG_LEN K_MIN K_MAX K_STEP PRESET
export RED GREEN YELLOW BLUE CYAN MAGENTA NC

# ============================================================================
# SCRIPT PRINCIPAL
# ============================================================================

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  Ensamblaje Metagenómico - MEGAHIT${NC}"
echo -e "${GREEN}  (Optimizado con GNU Parallel)${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo -e "Directorio de entrada: ${YELLOW}${INPUT_DIR}${NC}"
echo -e "Directorio de salida: ${YELLOW}${OUTPUT_DIR}${NC}"
echo -e "Threads por muestra: ${YELLOW}${THREADS_PER_SAMPLE}${NC}"
echo -e "Muestras en paralelo: ${YELLOW}${MAX_PARALLEL_JOBS}${NC}"
echo -e "K-mer range: ${YELLOW}${K_MIN}-${K_MAX} (step ${K_STEP})${NC}"
echo -e "Min contig length: ${YELLOW}${MIN_CONTIG_LEN} bp${NC}"
if [ -n "${PRESET}" ]; then
    echo -e "Preset: ${YELLOW}${PRESET}${NC}"
fi
echo ""

# ============================================================================
# VERIFICACIONES PREVIAS
# ============================================================================

# Verificar directorio de entrada
if [ ! -d "${INPUT_DIR}" ]; then
    echo -e "${RED}Error: El directorio de entrada no existe: ${INPUT_DIR}${NC}"
    exit 1
fi

# Verificar que MEGAHIT está instalado
if ! command -v megahit &> /dev/null; then
    echo -e "${RED}Error: MEGAHIT no está instalado${NC}"
    echo "Instalar con: conda install -c bioconda megahit"
    exit 1
fi

# Verificar versión de MEGAHIT
MEGAHIT_VERSION=$(megahit --version 2>&1 | grep -oP 'v\d+\.\d+\.\d+' | head -1)
echo -e "${CYAN}Versión de MEGAHIT: ${MEGAHIT_VERSION}${NC}"
echo ""

# Verificar que Python3 está disponible
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}Error: Python3 no está instalado${NC}"
    exit 1
fi

# Verificar que GNU Parallel está instalado
if ! command -v parallel &> /dev/null; then
    echo -e "${YELLOW}Advertencia: GNU Parallel no está instalado${NC}"
    echo -e "${YELLOW}Se procesarán las muestras secuencialmente${NC}"
    echo -e "${YELLOW}Para instalar: sudo apt-get install parallel${NC}"
    USE_PARALLEL=false
else
    USE_PARALLEL=true
fi

# Crear directorio de salida
mkdir -p "${OUTPUT_DIR}"

# Advertencia sobre recursos
echo -e "${MAGENTA}========================================${NC}"
echo -e "${MAGENTA}  ADVERTENCIA: RECURSOS${NC}"
echo -e "${MAGENTA}========================================${NC}"
echo -e "${YELLOW}MEGAHIT requiere recursos significativos:${NC}"
echo -e "${YELLOW}  - RAM: ~10-50 GB por muestra (depende del tamaño)${NC}"
echo -e "${YELLOW}  - CPU: Se beneficia de muchos threads (16+ recomendado)${NC}"
echo -e "${YELLOW}  - Tiempo: 1-12 horas por muestra (depende de complejidad)${NC}"
echo -e "${YELLOW}  - Disco: ~2-5× el tamaño de los FASTQ de entrada${NC}"
echo ""
echo -e "${YELLOW}Configuración actual:${NC}"
echo -e "${YELLOW}  - ${MAX_PARALLEL_JOBS} muestra(s) en paralelo${NC}"
echo -e "${YELLOW}  - ${THREADS_PER_SAMPLE} threads por muestra${NC}"
echo -e "${YELLOW}  - Uso total estimado: $((MAX_PARALLEL_JOBS * THREADS_PER_SAMPLE)) threads${NC}"
echo ""
echo -e "${CYAN}Verifica recursos disponibles:${NC}"
free -h | grep -E "Mem|Swap"
echo ""
df -h "${OUTPUT_DIR}" | tail -1
echo ""
echo -e "${MAGENTA}========================================${NC}"
echo ""

read -p "¿Continuar con el ensamblaje? (y/n): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo -e "${YELLOW}Ensamblaje cancelado por el usuario${NC}"
    exit 0
fi
echo ""

# ============================================================================
# PREPARAR LISTA DE MUESTRAS
# ============================================================================

# Crear archivo temporal con lista de muestras
SAMPLE_LIST=$(mktemp)

# Buscar archivos en subdirectorios de cada muestra
for SAMPLE_DIR in "${INPUT_DIR}"/*/; do
    if [ -d "${SAMPLE_DIR}" ]; then
        SAMPLE=$(basename "${SAMPLE_DIR}")
        R1_FILE="${SAMPLE_DIR}/${SAMPLE}${SUFFIX_R1}"
        R2_FILE="${SAMPLE_DIR}/${SAMPLE}${SUFFIX_R2}"
        
        # Verificar que existen ambos archivos
        if [ -f "${R1_FILE}" ] && [ -f "${R2_FILE}" ]; then
            # Agregar a la lista: R1_path|SAMPLE|R2_path
            echo "${R1_FILE}|${SAMPLE}|${R2_FILE}" >> "${SAMPLE_LIST}"
        else
            echo -e "${YELLOW}Advertencia: Archivos incompletos para ${SAMPLE}${NC}"
            [ ! -f "${R1_FILE}" ] && echo -e "${YELLOW}  Falta: ${R1_FILE}${NC}"
            [ ! -f "${R2_FILE}" ] && echo -e "${YELLOW}  Falta: ${R2_FILE}${NC}"
        fi
    fi
done

# Contar muestras
TOTAL_SAMPLES=$(wc -l < "${SAMPLE_LIST}")

if [ ${TOTAL_SAMPLES} -eq 0 ]; then
    echo -e "${RED}Error: No se encontraron muestras en ${INPUT_DIR}${NC}"
    echo -e "${RED}Buscando archivos con sufijos: ${SUFFIX_R1} y ${SUFFIX_R2}${NC}"
    rm -f "${SAMPLE_LIST}"
    exit 1
fi

echo -e "${GREEN}Se ensamblarán ${TOTAL_SAMPLES} muestras${NC}"
echo ""

# ============================================================================
# PROCESAMIENTO
# ============================================================================

GLOBAL_START_TIME=$(date +%s)

if [ "${USE_PARALLEL}" = true ] && [ ${MAX_PARALLEL_JOBS} -gt 1 ]; then
    # Procesamiento en paralelo con GNU Parallel
    echo -e "${GREEN}Iniciando ensamblaje en paralelo...${NC}"
    echo ""
    
    cat "${SAMPLE_LIST}" | parallel \
        --colsep '|' \
        --jobs ${MAX_PARALLEL_JOBS} \
        --progress \
        --joblog "${OUTPUT_DIR}/parallel_joblog.txt" \
        process_megahit {1} {2} {3} ${THREADS_PER_SAMPLE}
    
    PARALLEL_EXIT=$?
else
    # Procesamiento secuencial
    echo -e "${CYAN}Procesamiento secuencial...${NC}"
    echo ""
    
    COUNTER=0
    while IFS='|' read -r R1_FILE SAMPLE R2_FILE; do
        ((COUNTER++))
        echo -e "${MAGENTA}========================================${NC}"
        echo -e "${MAGENTA}  MUESTRA ${COUNTER} de ${TOTAL_SAMPLES}${NC}"
        echo -e "${MAGENTA}========================================${NC}"
        process_megahit "${R1_FILE}" "${SAMPLE}" "${R2_FILE}" ${THREADS_PER_SAMPLE}
        echo ""
    done < "${SAMPLE_LIST}"
    
    PARALLEL_EXIT=0
fi

GLOBAL_END_TIME=$(date +%s)
GLOBAL_ELAPSED=$((GLOBAL_END_TIME - GLOBAL_START_TIME))
GLOBAL_HOURS=$((GLOBAL_ELAPSED / 3600))
GLOBAL_MIN=$(((GLOBAL_ELAPSED % 3600) / 60))
GLOBAL_SEC=$((GLOBAL_ELAPSED % 60))

# ============================================================================
# RESUMEN FINAL CONSOLIDADO
# ============================================================================

echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  RESUMEN FINAL - TODOS LOS ENSAMBLAJES${NC}"
echo -e "${GREEN}========================================${NC}"
echo -e "Total de muestras procesadas: ${YELLOW}${TOTAL_SAMPLES}${NC}"
echo -e "Tiempo total: ${YELLOW}${GLOBAL_HOURS}h ${GLOBAL_MIN}m ${GLOBAL_SEC}s${NC}"
echo ""
echo -e "${GREEN}Los ensamblajes se encuentran en:${NC}"
echo -e "${YELLOW}${OUTPUT_DIR}${NC}"
echo ""

# Generar resumen consolidado
SUMMARY_FILE="${OUTPUT_DIR}/00_RESUMEN_GENERAL.txt"
echo "============================================================" > "${SUMMARY_FILE}"
echo "RESUMEN GENERAL - ENSAMBLAJES MEGAHIT" >> "${SUMMARY_FILE}"
echo "============================================================" >> "${SUMMARY_FILE}"
echo "Fecha: $(date +%Y-%m-%d\ %H:%M:%S)" >> "${SUMMARY_FILE}"
echo "Total de muestras: ${TOTAL_SAMPLES}" >> "${SUMMARY_FILE}"
echo "Tiempo total de procesamiento: ${GLOBAL_HOURS}h ${GLOBAL_MIN}m ${GLOBAL_SEC}s" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"
echo "PARÁMETROS:" >> "${SUMMARY_FILE}"
echo "  K-mer range: ${K_MIN}-${K_MAX} (step ${K_STEP})" >> "${SUMMARY_FILE}"
echo "  Min contig length: ${MIN_CONTIG_LEN} bp" >> "${SUMMARY_FILE}"
echo "  Preset: ${PRESET}" >> "${SUMMARY_FILE}"
echo "  Threads per sample: ${THREADS_PER_SAMPLE}" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"
echo "RESULTADOS POR MUESTRA:" >> "${SUMMARY_FILE}"
echo "============================================================" >> "${SUMMARY_FILE}"
printf "%-20s %12s %12s %12s %12s\n" "MUESTRA" "CONTIGS" "TOTAL (Mbp)" "N50 (bp)" "LONGEST (bp)" >> "${SUMMARY_FILE}"
echo "------------------------------------------------------------" >> "${SUMMARY_FILE}"

for SAMPLE_DIR in "${OUTPUT_DIR}"/*/; do
    if [ -d "${SAMPLE_DIR}" ] && [ -f "${SAMPLE_DIR}/assembly_stats.txt" ]; then
        SAMPLE=$(basename "${SAMPLE_DIR}")
        NUM_CONTIGS=$(grep "Número de contigs:" "${SAMPLE_DIR}/assembly_stats.txt" | grep -oP '\d+' | head -1 | sed 's/,//g')
        TOTAL_LENGTH=$(grep "Longitud total:" "${SAMPLE_DIR}/assembly_stats.txt" | grep -oP '\d+\.\d+(?= Mbp)' | head -1)
        N50=$(grep "N50:" "${SAMPLE_DIR}/assembly_stats.txt" | grep -oP '\d+' | head -1 | sed 's/,//g')
        LONGEST=$(grep "Contig más largo:" "${SAMPLE_DIR}/assembly_stats.txt" | grep -oP '\d+' | head -1 | sed 's/,//g')
        
        printf "%-20s %12s %12s %12s %12s\n" "${SAMPLE}" "${NUM_CONTIGS}" "${TOTAL_LENGTH}" "${N50}" "${LONGEST}" >> "${SUMMARY_FILE}"
    fi
done

echo "" >> "${SUMMARY_FILE}"
echo "============================================================" >> "${SUMMARY_FILE}"
echo "ARCHIVOS GENERADOS POR MUESTRA:" >> "${SUMMARY_FILE}"
echo "  - final.contigs.fa (contigs ensamblados)" >> "${SUMMARY_FILE}"
echo "  - assembly_stats.txt (estadísticas detalladas)" >> "${SUMMARY_FILE}"
echo "  - log (archivo de log de MEGAHIT)" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"
echo "============================================================" >> "${SUMMARY_FILE}"

echo -e "${GREEN}Resumen guardado en: ${YELLOW}${SUMMARY_FILE}${NC}"
echo ""

# Mostrar tabla en consola
echo -e "${CYAN}Tabla de resultados:${NC}"
column -t "${SUMMARY_FILE}" | tail -n +9 | head -n $((TOTAL_SAMPLES + 3))
echo ""

echo -e "${GREEN}Proceso completado - $(date +%Y-%m-%d\ %H:%M:%S)${NC}"

# Limpiar archivo temporal
rm -f "${SAMPLE_LIST}"

exit ${PARALLEL_EXIT}

