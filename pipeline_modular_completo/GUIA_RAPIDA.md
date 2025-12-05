# Guía Rápida: Pipeline Metagenómico Completo

## 🚀 Inicio Rápido

### 1. Preparar el Entorno

```bash
# Navegar al directorio del proyecto
cd /data2/shaday/prueba/

# Verificar que el pipeline esté disponible
ls pipeline_modular/metagenomics_pipeline.sh
```

### 2. Instalar Ambientes (Primera vez)

```bash
# Ejecutar script de instalación
bash pipeline_modular/INSTALACION_AMBIENTES.sh
```

Este script instalará los ambientes de micromamba para los módulos 6-9:
- `gtdbtk` (Taxonomía de bins)
- `prokka` (Anotación)
- `rgi` (Resistencia a antibióticos)
- `antismash` (Metabolitos secundarios)

### 3. Configurar Variables de Entorno

```bash
# Configurar GTDB-Tk (requerido para módulo 6)
export GTDBTK_DATA_PATH="/home_local/camda/gtdbtk_db"

# Opcional: Agregar al .bashrc para persistencia
echo 'export GTDBTK_DATA_PATH="/home_local/camda/gtdbtk_db"' >> ~/.bashrc
```

### 4. Ejecutar el Pipeline

```bash
# Ejecutar el pipeline interactivo
bash pipeline_modular/metagenomics_pipeline.sh
```

---

## 📋 Menú Interactivo

Al ejecutar el pipeline, verás el siguiente menú:

```
============================================================================
  PIPELINE METAGENÓMICO MODULAR
============================================================================

Seleccione los módulos a ejecutar:

  1. QC & Trimming (Trim Galore)
  2. Host Removal (Bowtie2)
  3. Assembly (MEGAHIT)
  4. Binning (MetaBAT2/MaxBin2/CONCOCT)
  5. Taxonomía de Reads (Kraken2)
  6. Taxonomía de Bins (GTDB-Tk)          ← NUEVO
  7. Anotación (Prokka)                   ← NUEVO
  8. Resistencia a Antibióticos (RGI)     ← NUEVO
  9. Metabolitos Secundarios (AntiSMASH)  ← NUEVO
  10. Análisis y Reportes

  A. Seleccionar TODOS los módulos
  R. Revisar selección actual
  E. Ejecutar pipeline con módulos seleccionados
  Q. Salir
```

---

## 🎯 Casos de Uso Comunes

### Caso 1: Ejecutar Pipeline Completo

**Escenario:** Primera vez analizando muestras, quieres ejecutar todos los módulos.

**Pasos:**
1. Presiona `A` (Seleccionar todos)
2. Presiona `E` (Ejecutar)
3. Selecciona modo Kraken2:
   - `1` para Simple (solo GTDB)
   - `2` para Dual (GTDB + PlusPFP) ← **Recomendado**
   - `3` para Triple (máxima cobertura)
4. Espera a que termine (puede tardar horas/días)

### Caso 2: Solo Análisis de Bins (Módulos 6-9)

**Escenario:** Ya tienes bins generados (módulo 4 completado), solo quieres analizarlos.

**Pasos:**
1. Presiona `6` (GTDB-Tk)
2. Presiona `7` (Prokka)
3. Presiona `8` (RGI)
4. Presiona `9` (AntiSMASH)
5. Presiona `R` (Revisar selección)
6. Presiona `E` (Ejecutar)

### Caso 3: Solo Taxonomía y Anotación

**Escenario:** Quieres clasificar y anotar bins, pero no analizar resistencia ni metabolitos.

**Pasos:**
1. Presiona `6` (GTDB-Tk)
2. Presiona `7` (Prokka)
3. Presiona `E` (Ejecutar)

### Caso 4: Agregar Análisis de Resistencia a Bins Existentes

**Escenario:** Ya ejecutaste Prokka (módulo 7), ahora quieres analizar resistencia.

**Pasos:**
1. Presiona `8` (RGI)
2. Presiona `E` (Ejecutar)

---

## 🔍 Verificar Resultados

### Estructura de Directorios

```bash
# Ver estructura de salida
tree -L 2 output/
```

**Salida esperada:**
```
output/
├── trim/                    # Módulo 1
├── host_removed/            # Módulo 2
├── megahit_assemblies/      # Módulo 3
├── binning/                 # Módulo 4
├── kraken2_dual/            # Módulo 5
├── gtdbtk/                  # Módulo 6 ✨
├── prokka/                  # Módulo 7 ✨
├── rgi/                     # Módulo 8 ✨
├── antismash/               # Módulo 9 ✨
└── analysis/                # Módulo 10
```

### Verificar Resultados por Módulo

#### Módulo 6: GTDB-Tk

```bash
# Ver clasificaciones de bacterias
cat output/gtdbtk/SRR5936076/gtdbtk.bac120.summary.tsv

# Ver clasificaciones de arqueas
cat output/gtdbtk/SRR5936076/gtdbtk.ar53.summary.tsv

# Ver log
tail output/gtdbtk/SRR5936076/gtdbtk.log
```

#### Módulo 7: Prokka

```bash
# Listar bins anotados
ls output/prokka/SRR5936076/

# Ver estadísticas de un bin
cat output/prokka/SRR5936076/bin.1/bin.1.txt

# Ver genes anotados
head output/prokka/SRR5936076/bin.1/bin.1.tsv
```

#### Módulo 8: RGI

```bash
# Ver genes de resistencia detectados
cat output/rgi/SRR5936076/bin.1/bin.1.txt

# Contar genes de resistencia por muestra
for sample in output/rgi/*/; do
    echo "=== $(basename $sample) ==="
    find $sample -name "*.txt" -exec wc -l {} \; | awk '{sum+=$1} END {print "Total genes:", sum-NR}'
done
```

#### Módulo 9: AntiSMASH

```bash
# Abrir reporte HTML de un bin
firefox output/antismash/SRR5936076/bin.1/index.html

# Listar clusters detectados
grep "Region" output/antismash/SRR5936076/bin.1/index.html
```

---

## ⚙️ Configuración Avanzada

### Ajustar Recursos Computacionales

Edita el archivo `metagenomics_pipeline.sh`:

```bash
# Líneas 56-61
THREADS_QC=20              # Hilos para QC/Trimming
THREADS_HOST=12            # Hilos para remoción de host
THREADS_ASSEMBLY=60        # Hilos para ensamblaje
THREADS_BINNING=40         # Hilos para binning y GTDB-Tk
THREADS_KRAKEN=40          # Hilos para Kraken2
```

### Cambiar Directorios

```bash
# Líneas 36-38
PROJECT_DIR="/files/shaday/4_cienegas"        # Directorio del proyecto
INPUT_DIR="${PROJECT_DIR}/MergedFastq"        # Directorio de entrada
OUTPUT_DIR="${PROJECT_DIR}/output"            # Directorio de salida
```

### Configurar Bases de Datos

```bash
# Líneas 49-52
BOWTIE2_INDEX="/path/to/bowtie2/index"        # Índice Bowtie2
KRAKEN2_GTDB="/path/to/k2_gtdb"               # Base de datos GTDB
KRAKEN2_PLUSPFP="/path/to/k2_pluspfp"         # Base de datos PlusPFP
KRAKEN2_EUPATH="/path/to/k2_eupathdb"         # Base de datos EuPathDB
```

---

## 🐛 Solución de Problemas

### Problema: "GTDB-Tk no está disponible"

**Solución:**
```bash
# Verificar que el ambiente existe
micromamba env list | grep gtdbtk

# Activar manualmente y probar
micromamba activate gtdbtk
gtdbtk --version
micromamba deactivate
```

### Problema: "Variable GTDBTK_DATA_PATH no configurada"

**Solución:**
```bash
# Configurar la variable
export GTDBTK_DATA_PATH="/home_local/camda/gtdbtk_db"

# Verificar que el directorio existe
ls -lh $GTDBTK_DATA_PATH
```

### Problema: "Base de datos CARD no cargada"

**Solución:**
```bash
# Activar ambiente RGI
micromamba activate rgi

# Descargar y cargar CARD
wget https://card.mcmaster.ca/latest/data -O card_data.tar.bz2
tar -xjf card_data.tar.bz2
rgi load --card_json card.json

# Verificar
rgi database --version

micromamba deactivate
```

### Problema: "No se encontró el directorio de binning"

**Solución:**
```bash
# Verificar que el módulo 4 se ejecutó correctamente
ls -lh output/binning/

# Si no existe, ejecutar módulo 4 primero
# Desde el menú: seleccionar opción 4 y ejecutar
```

### Problema: "Error en muestra ${SAMPLE}"

**Solución:**
```bash
# Ver el log específico del error
tail -n 50 output/gtdbtk/${SAMPLE}/gtdbtk.log
# o
tail -n 50 output/prokka/${SAMPLE}/bin.1/prokka.log

# Buscar mensajes de error
grep -i "error\|fail\|exception" output/gtdbtk/${SAMPLE}/gtdbtk.log
```

---

## 📊 Interpretación de Resultados

### GTDB-Tk (Módulo 6)

**Archivo:** `gtdbtk.bac120.summary.tsv`

**Columnas importantes:**
- `user_genome`: Nombre del bin
- `classification`: Taxonomía completa (d__Bacteria;p__...;g__...;s__...)
- `fastani_reference`: Genoma de referencia más cercano
- `fastani_ani`: Identidad ANI (>95% = misma especie)

**Ejemplo:**
```
bin.1   d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__Pseudomonas_aeruginosa
```

### Prokka (Módulo 7)

**Archivos importantes:**
- `.faa`: Secuencias de proteínas (para RGI)
- `.gbk`: Archivo GenBank (para AntiSMASH)
- `.gff`: Anotaciones (para visualización)
- `.tsv`: Tabla de genes

**Estadísticas en `.txt`:**
```
organism: Bacteria
contigs: 145
bases: 4523891
CDS: 4234
rRNA: 3
tRNA: 45
```

### RGI (Módulo 8)

**Archivo:** `bin.1.txt`

**Columnas importantes:**
- `Best_Hit_ARO`: Gen de resistencia detectado
- `Drug Class`: Clase de antibiótico
- `Resistance Mechanism`: Mecanismo de resistencia
- `% Identity`: Identidad con gen de referencia

**Ejemplo:**
```
mecA    beta-lactam    target alteration    98.5%
```

### AntiSMASH (Módulo 9)

**Archivo:** `index.html` (abrir en navegador)

**Información mostrada:**
- **Regiones detectadas:** Clusters biosintéticos
- **Tipo de cluster:** NRPS, PKS, terpene, bacteriocin, etc.
- **Genes involucrados:** Enzimas biosintéticas
- **Comparación:** Similitud con clusters conocidos

---

## 📈 Tiempos de Ejecución Estimados

**Hardware de referencia:** 40-60 cores, 256 GB RAM

| Módulo | Herramienta | Tiempo (2 muestras) | Tiempo (10 muestras) |
|--------|-------------|---------------------|----------------------|
| 1      | Trim Galore | 30 min              | 2-3 horas            |
| 2      | Bowtie2     | 1 hora              | 5 horas              |
| 3      | MEGAHIT     | 2-4 horas           | 10-20 horas          |
| 4      | MetaBAT2    | 1 hora              | 5 horas              |
| 5      | Kraken2     | 30 min              | 2 horas              |
| **6**  | **GTDB-Tk** | **4-8 horas**       | **20-40 horas**      |
| **7**  | **Prokka**  | **2-4 horas**       | **10-20 horas**      |
| **8**  | **RGI**     | **1-2 horas**       | **5-10 horas**       |
| **9**  | **AntiSMASH** | **4-8 horas**     | **20-40 horas**      |
| 10     | Análisis    | 10 min              | 30 min               |

**Total estimado:** 15-30 horas para 2 muestras, 80-150 horas para 10 muestras

---

## 🎓 Recursos Adicionales

### Documentación de Herramientas

- **GTDB-Tk:** https://ecogenomics.github.io/GTDBTk/
- **Prokka:** https://github.com/tseemann/prokka
- **RGI:** https://github.com/arpcard/rgi
- **AntiSMASH:** https://docs.antismash.secondarymetabolites.org/

### Bases de Datos

- **GTDB:** https://gtdb.ecogenomic.org/
- **CARD:** https://card.mcmaster.ca/
- **MIBiG:** https://mibig.secondarymetabolites.org/

### Tutoriales

- **GTDB-Tk Tutorial:** https://ecogenomics.github.io/GTDBTk/tutorials/
- **Prokka Tutorial:** https://github.com/tseemann/prokka#quick-start
- **RGI Tutorial:** https://github.com/arpcard/rgi#usage
- **AntiSMASH Tutorial:** https://docs.antismash.secondarymetabolites.org/understanding_output/

---

## ✅ Checklist de Verificación

Antes de ejecutar el pipeline completo:

- [ ] Micromamba instalado y en PATH
- [ ] Ambientes creados (gtdbtk, prokka, rgi, antismash)
- [ ] Variable GTDBTK_DATA_PATH configurada
- [ ] Base de datos GTDB descargada (~85 GB)
- [ ] Base de datos CARD cargada en RGI
- [ ] Datos de entrada en formato correcto (FASTQ pareados)
- [ ] Suficiente espacio en disco (>500 GB recomendado)
- [ ] Suficiente RAM (>128 GB recomendado)
- [ ] Permisos de escritura en directorio de salida

---

## 📞 Contacto

Para preguntas o problemas:
1. Revisar logs en `output/${modulo}/${muestra}/`
2. Consultar documentación oficial de cada herramienta
3. Verificar configuración de ambientes y bases de datos

---

**Última actualización:** 29 de noviembre de 2025  
**Versión del pipeline:** 1.0 (10 módulos completos)
