# Resumen de Implementación: Módulos 6-9 del Pipeline Metagenómico

## Estado Actual

Se han implementado exitosamente los **4 módulos faltantes** (6, 7, 8 y 9) en el pipeline metagenómico modular. Ahora el pipeline está completo con los 10 módulos funcionales.

---

## Módulos Implementados

### **Módulo 6: GTDB-Tk (Taxonomía de Bins)**

**Función:** `module_06_gtdbtk()`

**Descripción:** Clasifica taxonómicamente los bins generados por MetaBAT2 usando GTDB-Tk (Genome Taxonomy Database Toolkit).

**Características:**
- Activa el ambiente `gtdbtk` de micromamba
- Verifica la existencia del directorio de binning (módulo 4)
- Comprueba que GTDB-Tk esté instalado y la base de datos configurada (`GTDBTK_DATA_PATH`)
- Procesa cada muestra individualmente
- Cuenta los bins disponibles antes de procesar
- Ejecuta `gtdbtk classify_wf` con:
  - `--skip_ani_screen`: Desactiva el screening ANI (más rápido)
  - `--extension fa`: Extensión de archivos de bins
  - `--cpus ${THREADS_BINNING}`: Usa los hilos configurados para binning
- Genera logs individuales por muestra
- Produce archivos de salida:
  - `gtdbtk.bac120.summary.tsv`: Clasificación de bacterias
  - `gtdbtk.ar53.summary.tsv` o `gtdbtk.ar122.summary.tsv`: Clasificación de arqueas

**Dependencias:**
- Requiere Módulo 4 (Binning) completado
- Variable de entorno: `GTDBTK_DATA_PATH`
- Base de datos: GTDB release (r207 o superior)

**Salida:** `${OUTPUT_DIR}/gtdbtk/${SAMPLE}/`

---

### **Módulo 7: Prokka (Anotación Funcional)**

**Función:** `module_07_prokka()`

**Descripción:** Anota funcionalmente los bins usando Prokka, identificando genes, proteínas y características genómicas.

**Características:**
- Activa el ambiente `prokka` de micromamba
- Verifica directorios de binning y GTDB-Tk
- Si GTDB-Tk no está disponible, usa "Bacteria" como reino por defecto
- Determina el reino (Bacteria/Archaea) de cada bin basándose en resultados de GTDB-Tk
- Procesa cada bin individualmente con Prokka:
  - `--metagenome`: Modo metagenómico
  - `--kingdom`: Bacteria o Archaea según clasificación
  - `--mincontiglen 200`: Longitud mínima de contig
  - `--cpus 4`: 4 hilos por bin (ajustable)
  - `--force`: Sobrescribe resultados previos
- Genera múltiples archivos de salida por bin:
  - `.faa`: Secuencias de proteínas
  - `.ffn`: Secuencias de genes
  - `.gbk`: Archivo GenBank
  - `.gff`: Anotaciones en formato GFF
  - `.tsv`: Tabla de características

**Dependencias:**
- Requiere Módulo 4 (Binning) completado
- Opcional: Módulo 6 (GTDB-Tk) para clasificación precisa

**Salida:** `${OUTPUT_DIR}/prokka/${SAMPLE}/${bin_name}/`

---

### **Módulo 8: RGI (Resistencia a Antibióticos)**

**Función:** `module_08_rgi()`

**Descripción:** Analiza genes de resistencia a antibióticos usando RGI (Resistance Gene Identifier) y la base de datos CARD.

**Características:**
- Activa el ambiente `rgi` de micromamba
- Verifica que Prokka haya completado (necesita archivos `.faa`)
- Comprueba que la base de datos CARD esté cargada
- Procesa archivos de proteínas (`.faa`) de cada bin
- Ejecuta `rgi main` con:
  - `--input_type protein`: Entrada de proteínas
  - `--num_threads 4`: 4 hilos por bin
  - `--clean`: Limpia archivos temporales
- Genera archivos de resultados:
  - `.txt`: Tabla de genes de resistencia detectados
  - `.json`: Resultados en formato JSON

**Dependencias:**
- Requiere Módulo 7 (Prokka) completado
- Base de datos: CARD (Comprehensive Antibiotic Resistance Database)
- Comando para descargar CARD: `rgi load --card_json /path/to/card.json`

**Salida:** `${OUTPUT_DIR}/rgi/${SAMPLE}/${bin_name}/`

---

### **Módulo 9: AntiSMASH (Metabolitos Secundarios)**

**Función:** `module_09_antismash()`

**Descripción:** Identifica clusters de genes biosintéticos de metabolitos secundarios usando AntiSMASH.

**Características:**
- Activa el ambiente `antismash` de micromamba
- Verifica que Prokka haya completado (necesita archivos `.gbk`)
- Procesa archivos GenBank (`.gbk`) de cada bin
- Ejecuta `antismash` con:
  - `--genefinding-tool prodigal`: Usa Prodigal para predicción de genes
  - `--taxon bacteria`: Especifica taxón bacteriano
  - `--minlength 1000`: Longitud mínima de región
  - `--cpus 4`: 4 hilos por bin
- Genera resultados HTML interactivos con:
  - Clusters biosintéticos detectados
  - Genes involucrados
  - Comparación con clusters conocidos
  - Visualizaciones de regiones genómicas

**Dependencias:**
- Requiere Módulo 7 (Prokka) completado
- Base de datos: AntiSMASH databases (se descargan automáticamente en la primera ejecución)

**Salida:** `${OUTPUT_DIR}/antismash/${SAMPLE}/${bin_name}/`

---

## Integración en el Pipeline

### Cambios Realizados

1. **Funciones agregadas** en `metagenomics_pipeline.sh`:
   - `module_06_gtdbtk()`
   - `module_07_prokka()`
   - `module_08_rgi()`
   - `module_09_antismash()`

2. **Ejecución secuencial** actualizada en `run_pipeline()`:
   ```bash
   [ "${MODULES_SELECTED[6]}" = "1" ] && module_06_gtdbtk
   [ "${MODULES_SELECTED[7]}" = "1" ] && module_07_prokka
   [ "${MODULES_SELECTED[8]}" = "1" ] && module_08_rgi
   [ "${MODULES_SELECTED[9]}" = "1" ] && module_09_antismash
   ```

3. **Menú principal** ya incluye las opciones 6-9:
   - Opción 6: Taxonomía de Bins (GTDB-Tk)
   - Opción 7: Anotación (Prokka)
   - Opción 8: Resistencia a Antibióticos (RGI)
   - Opción 9: Metabolitos Secundarios (AntiSMASH)

---

## Flujo de Trabajo Completo

```
1. QC & Trimming (Trim Galore)
   ↓
2. Host Removal (Bowtie2)
   ↓
3. Assembly (MEGAHIT)
   ↓
4. Binning (MetaBAT2)
   ↓
5. Taxonomía de Reads (Kraken2) ──┐
   ↓                                │
6. Taxonomía de Bins (GTDB-Tk)     │
   ↓                                │
7. Anotación (Prokka)               │
   ↓                                │
8. Resistencia (RGI)                │
   ↓                                │
9. Metabolitos (AntiSMASH)          │
   ↓                                │
10. Análisis y Reportes ←───────────┘
```

---

## Requisitos de Ambientes

### Ambientes de Micromamba Necesarios

Los siguientes ambientes deben estar creados:

1. **qc_assembly**: Módulos 1-3 (Trim Galore, Bowtie2, MEGAHIT, Samtools)
2. **binning**: Módulo 4 (MetaBAT2, Bowtie2)
3. **kraken2**: Módulo 5 (Kraken2)
4. **gtdbtk**: Módulo 6 (GTDB-Tk)
5. **prokka**: Módulo 7 (Prokka)
6. **rgi**: Módulo 8 (RGI)
7. **antismash**: Módulo 9 (AntiSMASH)
8. **analysis**: Módulo 10 (Python, biom-format, pandas, matplotlib)

### Comandos de Instalación

```bash
# Módulo 6: GTDB-Tk
micromamba create -n gtdbtk -c bioconda gtdbtk

# Módulo 7: Prokka
micromamba create -n prokka -c bioconda prokka

# Módulo 8: RGI
micromamba create -n rgi -c bioconda rgi

# Módulo 9: AntiSMASH
micromamba create -n antismash -c bioconda antismash
```

---

## Bases de Datos Requeridas

### GTDB-Tk (Módulo 6)
- **Ubicación:** Variable `GTDBTK_DATA_PATH`
- **Descarga:** https://data.gtdb.ecogenomic.org/
- **Tamaño:** ~85 GB (release r214)
- **Configuración:**
  ```bash
  export GTDBTK_DATA_PATH="/path/to/gtdbtk_db"
  ```

### CARD (Módulo 8)
- **Descarga:** https://card.mcmaster.ca/download
- **Instalación:**
  ```bash
  wget https://card.mcmaster.ca/latest/data
  tar -xvf data
  rgi load --card_json card.json
  ```

### AntiSMASH (Módulo 9)
- **Descarga automática:** Las bases de datos se descargan automáticamente en la primera ejecución
- **Ubicación:** `~/.local/share/antismash/`

---

## Verificación de Implementación

### Comprobar que los módulos están en el script

```bash
grep -n "^module_0[6-9]" metagenomics_pipeline.sh
```

**Salida esperada:**
```
598:module_06_gtdbtk() {
678:module_07_prokka() {
778:module_08_rgi() {
858:module_09_antismash() {
```

### Comprobar que se ejecutan en el pipeline

```bash
grep -n "MODULES_SELECTED\[6\]" metagenomics_pipeline.sh
grep -n "MODULES_SELECTED\[7\]" metagenomics_pipeline.sh
grep -n "MODULES_SELECTED\[8\]" metagenomics_pipeline.sh
grep -n "MODULES_SELECTED\[9\]" metagenomics_pipeline.sh
```

---

## Pruebas Recomendadas

### 1. Verificar sintaxis del script

```bash
bash -n metagenomics_pipeline.sh
```

### 2. Ejecutar solo el módulo 6 (GTDB-Tk)

Desde el menú interactivo:
1. Seleccionar opción `6`
2. Seleccionar opción `E` (Ejecutar)

### 3. Ejecutar módulos 6-9 en secuencia

Desde el menú interactivo:
1. Seleccionar opciones `6`, `7`, `8`, `9`
2. Seleccionar opción `E` (Ejecutar)

### 4. Ejecutar pipeline completo

Desde el menú interactivo:
1. Seleccionar opción `A` (Todos los módulos)
2. Configurar modo Kraken2
3. Seleccionar opción `E` (Ejecutar)

---

## Estructura de Salida

```
${OUTPUT_DIR}/
├── trim/                    # Módulo 1
├── host_removed/            # Módulo 2
├── megahit_assemblies/      # Módulo 3
├── binning/                 # Módulo 4
│   └── ${SAMPLE}/
│       └── metabat2/
│           ├── bin.1.fa
│           ├── bin.2.fa
│           └── ...
├── kraken2_${MODE}/         # Módulo 5
├── gtdbtk/                  # Módulo 6 ✨ NUEVO
│   └── ${SAMPLE}/
│       ├── gtdbtk.bac120.summary.tsv
│       ├── gtdbtk.ar53.summary.tsv
│       └── gtdbtk.log
├── prokka/                  # Módulo 7 ✨ NUEVO
│   └── ${SAMPLE}/
│       └── ${bin_name}/
│           ├── ${bin_name}.faa
│           ├── ${bin_name}.ffn
│           ├── ${bin_name}.gbk
│           ├── ${bin_name}.gff
│           └── prokka.log
├── rgi/                     # Módulo 8 ✨ NUEVO
│   └── ${SAMPLE}/
│       └── ${bin_name}/
│           ├── ${bin_name}.txt
│           ├── ${bin_name}.json
│           └── rgi.log
├── antismash/               # Módulo 9 ✨ NUEVO
│   └── ${SAMPLE}/
│       └── ${bin_name}/
│           ├── index.html
│           ├── regions.js
│           └── antismash.log
└── analysis/                # Módulo 10
    ├── combined_results.biom
    └── resultados_metagenomicos/
```

---

## Notas Importantes

### Recursos Computacionales

Los módulos 6-9 son **intensivos en recursos**:

- **GTDB-Tk:** Requiere ~150 GB de RAM y espacio temporal
- **Prokka:** ~2-4 GB RAM por bin
- **RGI:** ~1-2 GB RAM por bin
- **AntiSMASH:** ~4-8 GB RAM por bin, puede tardar horas

### Paralelización

Actualmente, los bins se procesan **secuencialmente** dentro de cada muestra. Para acelerar:

1. Aumentar `--cpus` en cada herramienta
2. Implementar procesamiento paralelo de bins (GNU Parallel)
3. Procesar múltiples muestras simultáneamente

### Manejo de Errores

Cada módulo:
- Verifica dependencias antes de ejecutar
- Registra errores en archivos `.log`
- Continúa con el siguiente bin/muestra si uno falla
- No detiene el pipeline completo por errores individuales

---

## Próximos Pasos

1. **Crear los ambientes de micromamba** para módulos 6-9
2. **Descargar las bases de datos** necesarias (GTDB, CARD)
3. **Probar cada módulo individualmente** con datos de prueba
4. **Ejecutar pipeline completo** en las muestras SRR5936076 y SRR5936077
5. **Optimizar recursos** según el hardware disponible
6. **Documentar resultados** y ajustar parámetros si es necesario

---

## Contacto y Soporte

Para problemas o preguntas sobre la implementación:
- Revisar logs en `${OUTPUT_DIR}/${module}/${SAMPLE}/`
- Verificar que los ambientes estén correctamente instalados
- Comprobar que las bases de datos estén configuradas

---

**Fecha de implementación:** 29 de noviembre de 2025  
**Versión del pipeline:** 1.0 (completo, 10 módulos)  
**Estado:** ✅ Implementación completa, listo para pruebas
