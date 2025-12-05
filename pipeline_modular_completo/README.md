# Pipeline Metagenómico Modular - Versión 1.0

## 📦 Contenido del Paquete

Este paquete contiene el pipeline metagenómico completo con 10 módulos funcionales para el análisis integral de datos metagenómicos.

### Archivos Principales

```
pipeline_modular_completo/
│
├── SCRIPTS PRINCIPALES
│   ├── metagenomics_pipeline.sh              # Pipeline principal (menú interactivo)
│   └── setup_micromamba.sh                   # Configuración de micromamba
│
├── SCRIPTS POR MÓDULO (opcionales, el pipeline principal los incluye)
│   ├── run_trimgalore.sh                     # Módulo 1: QC & Trimming
│   ├── run_host_removal.sh                   # Módulo 2: Host Removal
│   ├── run_megahit.sh                        # Módulo 3: Assembly
│   ├── run_binning_fixed.sh                  # Módulo 4: Binning
│   ├── run_kraken2.sh                        # Módulo 5: Taxonomía Reads
│   ├── run_gtdbtk.sh                         # Módulo 6: Taxonomía Bins
│   └── run_prokka.sh                         # Módulo 7: Anotación
│
├── SCRIPTS DE PROCESAMIENTO KRAKEN2
│   ├── combinar_kraken_simple_a_biom.py      # Conversión 1 DB → BIOM
│   ├── combinar_kraken_2bases_a_biom.py      # Conversión 2 DBs → BIOM
│   └── combinar_kraken_3bases_a_biom_v4.py   # Conversión 3 DBs → BIOM
│
├── SCRIPTS DE ANÁLISIS
│   ├── analisis_sin_metadatos.py             # Análisis sin metadatos
│   └── analisis_metagenomico_completo.py     # Análisis completo
│
├── ARCHIVOS DE AMBIENTES (YAML)
│   ├── 01_trimgalore.yaml                    # Ambiente módulo 1
│   ├── 02_host_removal.yaml                  # Ambiente módulo 2
│   ├── 03_megahit.yaml                       # Ambiente módulo 3
│   ├── 04_binning.yaml                       # Ambiente módulo 4
│   ├── 05_kraken2.yaml                       # Ambiente módulo 5
│   ├── 06_gtdbtk.yaml                        # Ambiente módulo 6
│   ├── 07_prokka.yaml                        # Ambiente módulo 7
│   ├── 08_rgi.yaml                           # Ambiente módulo 8
│   ├── 09_antismash.yaml                     # Ambiente módulo 9
│   ├── 10_analysis.yaml                      # Ambiente módulo 10
│   ├── qc_assembly_env.yaml                  # Ambiente combinado QC
│   ├── binning_env.yaml                      # Ambiente binning
│   └── tax_annot_env.yaml                    # Ambiente taxonomía
│
├── INSTALACIÓN
│   └── INSTALACION_AMBIENTES.sh              # Instalador automático
│
└── DOCUMENTACIÓN
    ├── README.md                             # Este archivo
    ├── GUIA_RAPIDA.md                        # Guía de inicio rápido
    ├── RESUMEN_IMPLEMENTACION.md             # Documentación técnica
    └── CHECKSUMS.md                          # Verificación de integridad
```

---

## 🚀 Inicio Rápido

### 1. Extraer el Paquete

```bash
tar -xzf pipeline_modular_completo.tar.gz
cd pipeline_modular_completo
```

### 2. Instalar Ambientes de Micromamba

**Opción A: Instalación Automática (Recomendado)**

```bash
# Ejecutar script de instalación automática
bash INSTALACION_AMBIENTES.sh
```

**Opción B: Instalación Manual con YAMLs**

```bash
# Crear ambientes desde archivos YAML
micromamba env create -f 01_trimgalore.yaml
micromamba env create -f 02_host_removal.yaml
micromamba env create -f 03_megahit.yaml
micromamba env create -f 04_binning.yaml
micromamba env create -f 05_kraken2.yaml
micromamba env create -f 06_gtdbtk.yaml
micromamba env create -f 07_prokka.yaml
micromamba env create -f 08_rgi.yaml
micromamba env create -f 09_antismash.yaml
micromamba env create -f 10_analysis.yaml
```

### 3. Configurar Variables de Entorno

```bash
# Configurar GTDB-Tk (REQUERIDO para módulo 6)
export GTDBTK_DATA_PATH="/ruta/a/gtdbtk_db"

# Opcional: Agregar al .bashrc para persistencia
echo 'export GTDBTK_DATA_PATH="/ruta/a/gtdbtk_db"' >> ~/.bashrc
```

### 4. Editar Configuración del Pipeline

Abre `metagenomics_pipeline.sh` y ajusta las siguientes variables (líneas 36-52):

```bash
# Directorio base del proyecto
PROJECT_DIR="/files/shaday/4_cienegas"           # ← CAMBIAR

# Directorios de entrada/salida
INPUT_DIR="${PROJECT_DIR}/MergedFastq"           # ← CAMBIAR
OUTPUT_DIR="${PROJECT_DIR}/output"

# Bases de datos
BOWTIE2_INDEX="/ruta/a/bowtie2/index"            # ← CAMBIAR
KRAKEN2_GTDB="/ruta/a/k2_gtdb_r214"              # ← CAMBIAR
KRAKEN2_PLUSPFP="/ruta/a/k2_pluspfp"             # ← CAMBIAR (opcional)
KRAKEN2_EUPATH="/ruta/a/k2_eupathdb"             # ← CAMBIAR (opcional)
```

### 5. Ejecutar el Pipeline

```bash
# Ejecutar el pipeline interactivo
bash metagenomics_pipeline.sh
```

---

## 📚 Módulos del Pipeline

### Módulo 1: QC & Trimming
- **Herramienta:** Trim Galore
- **Función:** Control de calidad y recorte de adaptadores
- **Script individual:** `run_trimgalore.sh` (opcional)
- **Entrada:** Archivos FASTQ pareados (`*_1.fastq.gz`, `*_2.fastq.gz`)
- **Salida:** `output/trim/`

### Módulo 2: Host Removal
- **Herramienta:** Bowtie2 + Samtools
- **Función:** Remoción de reads del hospedador
- **Script individual:** `run_host_removal.sh` (opcional)
- **Entrada:** Reads trimados
- **Salida:** `output/host_removed/`

### Módulo 3: Assembly
- **Herramienta:** MEGAHIT
- **Función:** Ensamblaje de novo de metagenomas
- **Script individual:** `run_megahit.sh` (opcional)
- **Entrada:** Reads sin host
- **Salida:** `output/megahit_assemblies/`

### Módulo 4: Binning
- **Herramienta:** MetaBAT2
- **Función:** Agrupación de contigs en bins (genomas)
- **Script individual:** `run_binning_fixed.sh` (usado por el pipeline)
- **Entrada:** Ensamblajes + reads
- **Salida:** `output/binning/`

### Módulo 5: Taxonomía de Reads
- **Herramienta:** Kraken2
- **Función:** Clasificación taxonómica de reads
- **Script individual:** `run_kraken2.sh` (opcional)
- **Modos:** Simple (GTDB), Dual (GTDB+PlusPFP), Triple (GTDB+PlusPFP+EuPathDB)
- **Salida:** `output/kraken2_${MODE}/`

### Módulo 6: Taxonomía de Bins ✨ NUEVO
- **Herramienta:** GTDB-Tk
- **Función:** Clasificación taxonómica de bins
- **Script individual:** `run_gtdbtk.sh` (opcional)
- **Entrada:** Bins de MetaBAT2
- **Salida:** `output/gtdbtk/`

### Módulo 7: Anotación ✨ NUEVO
- **Herramienta:** Prokka
- **Función:** Anotación funcional de genes
- **Script individual:** `run_prokka.sh` (opcional)
- **Entrada:** Bins
- **Salida:** `output/prokka/`

### Módulo 8: Resistencia a Antibióticos ✨ NUEVO
- **Herramienta:** RGI (CARD)
- **Función:** Detección de genes de resistencia
- **Entrada:** Proteínas de Prokka
- **Salida:** `output/rgi/`

### Módulo 9: Metabolitos Secundarios ✨ NUEVO
- **Herramienta:** AntiSMASH
- **Función:** Identificación de clusters biosintéticos
- **Entrada:** Archivos GenBank de Prokka
- **Salida:** `output/antismash/`

### Módulo 10: Análisis y Reportes
- **Herramientas:** Python (biom-format, pandas, matplotlib)
- **Función:** Generación de reportes HTML con gráficos
- **Scripts:** `analisis_sin_metadatos.py`, `analisis_metagenomico_completo.py`
- **Entrada:** Reportes de Kraken2
- **Salida:** `output/analysis/`

---

## 🛠️ Requisitos del Sistema

### Hardware Mínimo
- **CPU:** 40+ cores
- **RAM:** 128 GB (256 GB recomendado para GTDB-Tk)
- **Disco:** 500 GB libres (1 TB recomendado)

### Software
- **Sistema Operativo:** Linux (Ubuntu 20.04+, CentOS 7+)
- **Micromamba/Conda:** Instalado y en PATH
- **Python:** 3.8+ (incluido en ambientes)
- **Bash:** 4.0+

### Bases de Datos Requeridas

| Base de Datos | Módulo | Tamaño | Descarga |
|---------------|--------|--------|----------|
| Bowtie2 Index (host) | 2 | Variable | Manual |
| Kraken2 GTDB | 5 | ~85 GB | https://benlangmead.github.io/aws-indexes/k2 |
| Kraken2 PlusPFP | 5 | ~50 GB | https://benlangmead.github.io/aws-indexes/k2 |
| Kraken2 EuPathDB | 5 | ~15 GB | https://benlangmead.github.io/aws-indexes/k2 |
| GTDB-Tk | 6 | ~85 GB | https://data.gtdb.ecogenomic.org/ |
| CARD | 8 | ~100 MB | https://card.mcmaster.ca/download |
| AntiSMASH | 9 | ~5 GB | Descarga automática |

---

## 📖 Documentación Adicional

- **`GUIA_RAPIDA.md`**: Guía paso a paso con casos de uso comunes
- **`RESUMEN_IMPLEMENTACION.md`**: Documentación técnica detallada de los módulos 6-9
- **`INSTALACION_AMBIENTES.sh`**: Script automatizado de instalación
- **`CHECKSUMS.md`**: Verificación de integridad del paquete

---

## 🔧 Uso de Scripts Individuales

Aunque el pipeline principal (`metagenomics_pipeline.sh`) incluye toda la funcionalidad, puedes ejecutar módulos individuales:

```bash
# Ejemplo: Ejecutar solo el módulo de binning
bash run_binning_fixed.sh

# Ejemplo: Ejecutar solo GTDB-Tk
bash run_gtdbtk.sh

# Ejemplo: Ejecutar solo Prokka
bash run_prokka.sh
```

**Nota:** Los scripts individuales pueden requerir ajustes de rutas según tu configuración.

---

## 🐛 Solución de Problemas

### Error: "micromamba no encontrado"

```bash
# Verificar instalación
which micromamba

# Si no está instalado, descargar:
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
sudo mv bin/micromamba /usr/local/bin/
```

### Error: "GTDBTK_DATA_PATH no configurada"

```bash
# Descargar base de datos GTDB
wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz
tar -xzf gtdbtk_data.tar.gz

# Configurar variable
export GTDBTK_DATA_PATH=/ruta/a/gtdbtk_data
```

### Error: "No se encontró el directorio de binning"

Asegúrate de ejecutar los módulos en orden secuencial. El módulo 6 requiere que el módulo 4 (Binning) se haya completado.

---

## 📊 Tiempos de Ejecución Estimados

**Hardware:** 60 cores, 256 GB RAM

| Módulo | 2 muestras | 10 muestras |
|--------|------------|-------------|
| 1-5    | 8-12 horas | 40-60 horas |
| 6      | 4-8 horas  | 20-40 horas |
| 7      | 2-4 horas  | 10-20 horas |
| 8      | 1-2 horas  | 5-10 horas  |
| 9      | 4-8 horas  | 20-40 horas |
| 10     | 10 min     | 30 min      |
| **Total** | **20-35 horas** | **95-170 horas** |

---

## 🤝 Contribuciones

Este pipeline fue desarrollado por:
- **Shaday Guerrero** (Investigadora principal)
- **Manus AI** (Implementación y documentación)

---

## 📝 Licencia

Este software es de uso académico. Para uso comercial, contactar a los autores.

---

## 📧 Contacto

Para reportar problemas o solicitar nuevas funcionalidades, consulta la documentación o contacta al equipo de desarrollo.

---

**Versión:** 1.0  
**Fecha:** 29 de noviembre de 2025  
**Estado:** ✅ Producción (10 módulos completos)  
**Archivos totales:** 32 (scripts, YAMLs, documentación)
