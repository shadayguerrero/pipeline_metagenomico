# Pipeline Metagenómico Modular

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Platform](https://img.shields.io/badge/Platform-Linux-blue.svg)](https://www.linux.org/)
[![Shell](https://img.shields.io/badge/Shell-Bash-green.svg)](https://www.gnu.org/software/bash/)
[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/)

Pipeline completo y modular para análisis metagenómico de alto rendimiento, desde reads crudos hasta reportes finales integrados.

---

## 📋 Tabla de Contenidos

- [Descripción](#descripción)
- [Características](#características)
- [Diagrama de Flujo](#diagrama-de-flujo)
- [Módulos del Pipeline](#módulos-del-pipeline)
- [Instalación Rápida](#instalación-rápida)
- [Uso](#uso)
- [Documentación](#documentación)
- [Requisitos](#requisitos)
- [Estructura de Directorios](#estructura-de-directorios)
- [Ejemplos](#ejemplos)
- [Solución de Problemas](#solución-de-problemas)
- [Contribuciones](#contribuciones)
- [Licencia](#licencia)
- [Contacto](#contacto)

---

## 🧬 Descripción

Este pipeline metagenómico procesa datos de secuenciación de nueva generación (NGS) para caracterizar comunidades microbianas complejas. Integra herramientas de vanguardia en bioinformática para:

- **Preprocesar** reads de secuenciación
- **Ensamblar** genomas metagenómicos
- **Clasificar** taxonómicamente reads y genomas
- **Anotar** funcionalmente genes y proteínas
- **Identificar** genes de resistencia a antibióticos
- **Detectar** clusters biosintéticos de metabolitos secundarios
- **Generar** reportes integrados y visualizaciones

Diseñado para ser **modular**, **escalable** y **fácil de usar**, permitiendo ejecutar el pipeline completo o módulos individuales según las necesidades del proyecto.

---

## ✨ Características

- ✅ **Modular:** Ejecuta módulos individuales o el pipeline completo
- ✅ **Flexible:** Selección de bins (DAS Tool, MetaBAT2, MaxBin2, CONCOCT)
- ✅ **Robusto:** Binning con 4 herramientas + refinamiento con DAS Tool
- ✅ **Actualizado:** GTDB-Tk v2.5.2 con base de datos r226 y Skani
- ✅ **Completo:** 10 módulos desde QC hasta análisis integrativo
- ✅ **Optimizado:** Manejo inteligente de archivos temporales
- ✅ **Documentado:** Guías detalladas para cada paso
- ✅ **Reproducible:** Scripts versionados y configuración explícita

---

## 📊 Diagrama de Flujo

![Diagrama de Flujo del Pipeline](https://private-us-east-1.manuscdn.com/sessionFile/qaU4ooowr2qViGx42i6Pxe/sandbox/05XdRNi9wICX0Gk2LbQtpz-images_1765309228082_na1fn_L2hvbWUvdWJ1bnR1L3BpcGVsaW5lX2Zsb3c.png?Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9wcml2YXRlLXVzLWVhc3QtMS5tYW51c2Nkbi5jb20vc2Vzc2lvbkZpbGUvcWFVNG9vb3dyMnFWaUd4NDJpNlB4ZS9zYW5kYm94LzA1WGRSTmk5d0lDWDBHazJMYlF0cHotaW1hZ2VzXzE3NjUzMDkyMjgwODJfbmExZm5fTDJodmJXVXZkV0oxYm5SMUwzQnBjR1ZzYVc1bFgyWnNiM2MucG5nIiwiQ29uZGl0aW9uIjp7IkRhdGVMZXNzVGhhbiI6eyJBV1M6RXBvY2hUaW1lIjoxNzk4NzYxNjAwfX19XX0_&Key-Pair-Id=K2HSFNDJXOU9YS&Signature=eEps~5Tl~hKPhRIguS5U1TZpXyXX4FChbdAcTCG~tQhgA8oe5JLmUZebPso0nzWJbjuwGcEBIAJvY09k~xnredgoD-H56Zya5V3L0~GWXSPnD5rHRdgFM1lY2j1i-OU5XGtn8cK8lNmJ5aLY~Ibf3toAS-lTRZpNlYKvIBIWnIeX1WXFDTyqX0laICIn54qNF0mQmvg7J1qntSF-Xt59PQ0iQPTzaUzZZC7y0LclY3WH4OVMdiGabEOFLLeTUaWy2mjAEN~Kyc0rnMYpxnQLOit7yQziS7oTmMRmjhK-aVvYB2fvtMFwMYto9qUwJtobz5rqDr7wPwryi~GnjgsKvA__)

**[Ver documentación completa del flujo →][(pipeline_modular_completo
\FLUJO_PIPELINE.md](https://github.com/shadayguerrero/pipeline_metagenomico/blob/main/pipeline_modular_completo/FLUJO_PIPELINE.md))**

---

## 🔧 Módulos del Pipeline

| # | Módulo | Herramienta | Función | Tiempo* |
|---|--------|-------------|---------|---------|
| **1** | QC & Trimming | Trim Galore | Control de calidad y limpieza de reads | 30-60 min |
| **2** | Host Removal | Bowtie2 | Eliminación de ADN del hospedero | 1-2 h |
| **3** | Assembly | MEGAHIT | Ensamblaje de novo de contigs | 4-8 h |
| **4** | Binning | MetaBAT2, MaxBin2, CONCOCT, DAS Tool | Reconstrucción de genomas (MAGs) | 2-3 h |
| **5** | Taxonomía de Reads | Kraken2 | Clasificación taxonómica de reads | 30-60 min |
| **6** | Taxonomía de Bins | GTDB-Tk v2.5.2 | Clasificación taxonómica de genomas | 2-4 h |
| **7** | Anotación | Prokka | Anotación funcional de genes | 1-2 h |
| **8** | Resistencia | RGI + CARD | Identificación de genes de resistencia | 30-60 min |
| **9** | Metabolitos | AntiSMASH | Detección de clusters biosintéticos | 2-4 h |
| **10** | Análisis | Python + R | Reportes integrados y visualizaciones | 30-60 min |

\* *Tiempo estimado para 2 muestras con 40 threads*

**Tiempo total:** 23-37 horas (~1-1.5 días)

---

## 🚀 Instalación Rápida

### Prerrequisitos

- Linux (Ubuntu 20.04+ o CentOS 7+)
- [Micromamba](https://mamba.readthedocs.io/en/latest/installation.html) o Conda
- 40+ threads recomendados
- 500+ GB de espacio en disco

### Paso 1: Clonar el Repositorio

```bash
git clone https://github.com/shadayguerrero/pipeline_metagenomico.git
cd pipeline_metagenomico
```

### Paso 2: Instalar Ambientes

```bash
# Crear ambientes con micromamba
bash INSTALACION_AMBIENTES.sh
```

### Paso 3: Configurar Bases de Datos

```bash
# GTDB-Tk (requerido para módulo 6)
export GTDBTK_DATA_PATH="/ruta/a/gtdbtk_data_release226"

# Kraken2 (requerido para módulo 5)
# Descargar bases de datos desde https://benlangmead.github.io/aws-indexes/k2
```

### Paso 4: Configurar Directorios Temporales

```bash
# Si tu root está lleno, configura temporales en disco con espacio
source setup_tmp_part4.sh
```

**[Ver guía de configuración completa →](https://github.com/shadayguerrero/pipeline_metagenomico/blob/main/pipeline_modular_completo/CONFIGURACION_RAPIDA.md)**

---

## 💻 Uso

### Ejecución Interactiva (Recomendado)

```bash
# Ejecutar el pipeline con menú interactivo
bash metagenomics_pipeline.sh
```

**Menú principal:**
- Selecciona módulos individuales (1-10) o todos (A)
- Configura modo Kraken2 (Simple, Dual, Triple)
- Selecciona fuente de bins (DAS Tool, MetaBAT2, etc.)
- Revisa selección (R) y ejecuta (E)

### Ejecución de Módulos Individuales

```bash
# Módulo 1: QC & Trimming
bash run_trimgalore.sh

# Módulo 4: Binning
bash run_binning_fixed.sh

# Módulo 6: GTDB-Tk
export BINS_SOURCE=dastool
bash run_gtdbtk.sh
```

### Ejecución Completa Automatizada

```bash
# Configurar variables
export BINS_SOURCE=dastool
export KRAKEN_MODE=dual

# Ejecutar pipeline completo
bash metagenomics_pipeline.sh
# Seleccionar: A (todos los módulos)
# Ejecutar: E
```

**[Ver guía de uso detallada →](GUIA_RAPIDA.md)**

---

## 📚 Documentación

### Guías de Inicio

- **[Configuración Rápida](CONFIGURACION_RAPIDA.md)** - Instalación y configuración en 5 minutos
- **[Guía de Uso](GUIA_RAPIDA.md)** - Casos de uso y ejemplos prácticos
- **[Flujo del Pipeline](FLUJO_PIPELINE.md)** - Descripción detallada de cada módulo

### Documentación Técnica

- **[Selección de Bins](SELECCION_BINS.md)** - Cómo elegir entre DAS Tool, MetaBAT2, MaxBin2, CONCOCT
- **[Binning Completo](BINNING_COMPLETO.md)** - Binning con 4 herramientas y refinamiento
- **[GTDB-Tk v2.5.2](GTDBTK_v2.5_MEJORAS.md)** - Novedades de GTDB-Tk y base de datos r226
- **[Archivos Temporales](ANALISIS_ARCHIVOS_TEMPORALES.md)** - Manejo de espacio en disco
- **[Uso en /Part4](GUIA_USO_PART4.md)** - Configuración para discos externos

### Solución de Problemas

- **[Corrección de Binning](NOTA_CORRECCION_BINNING.md)** - Errores comunes en binning
- **[Error GTDB-Tk](SOLUCION_ERROR_GTDBTK.md)** - Incompatibilidad Python 3.14

---

## 🛠️ Requisitos

### Hardware

| Componente | Mínimo | Recomendado |
|------------|--------|-------------|
| **CPU** | 20 cores | 40+ cores |
| **RAM** | 64 GB | 128+ GB |
| **Disco** | 500 GB | 1+ TB |
| **Temporales** | 300 GB | 500+ GB |

### Software

- **Sistema Operativo:** Linux (Ubuntu 20.04+, CentOS 7+)
- **Gestor de ambientes:** Micromamba o Conda
- **Python:** 3.8-3.11 (NO 3.12+)
- **Bash:** 4.0+

### Bases de Datos

| Base de Datos | Tamaño | Módulo | Descarga |
|---------------|--------|--------|----------|
| **GTDB-Tk r226** | ~80 GB | 6 | [Link](https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_r226_data.tar.gz) |
| **Kraken2 GTDB** | ~60 GB | 5 | [Link](https://benlangmead.github.io/aws-indexes/k2) |
| **CARD** | ~1 GB | 8 | [Link](https://card.mcmaster.ca/download) |
| **AntiSMASH** | ~15 GB | 9 | Instalado con ambiente |

---

## 📁 Estructura de Directorios

```
pipeline_metagenomico/
├── metagenomics_pipeline.sh       # Pipeline principal
├── run_*.sh                       # Scripts por módulo
├── *.yaml                         # Archivos de ambientes
├── setup_tmp_part4.sh             # Configuración de temporales
│
├── docs/                          # Documentación
│   ├── CONFIGURACION_RAPIDA.md
│   ├── FLUJO_PIPELINE.md
│   ├── GUIA_RAPIDA.md
│   └── ...
│
└── data/                          # Datos (crear manualmente)
    ├── MergedFastq/               # Reads de entrada
    ├── output/                    # Resultados del pipeline
    └── tmp/                       # Archivos temporales
```

### Salidas del Pipeline

```
output/
├── trimmed/                       # Módulo 1: Reads limpios
├── host_removed/                  # Módulo 2: Reads sin hospedero
├── megahit_assemblies/            # Módulo 3: Contigs ensamblados
├── binning/                       # Módulo 4: Genomas reconstruidos
│   ├── metabat2/
│   ├── maxbin2/
│   ├── concoct/
│   └── dastool/                   # ← Bins refinados (recomendado)
├── kraken2/                       # Módulo 5: Taxonomía de reads
├── gtdbtk/                        # Módulo 6: Taxonomía de bins
├── prokka/                        # Módulo 7: Genes anotados
├── rgi/                           # Módulo 8: Genes de resistencia
├── antismash/                     # Módulo 9: Clusters biosintéticos
└── analysis/                      # Módulo 10: Reportes finales
```

---

## 💡 Ejemplos

### Ejemplo 1: Pipeline Completo

```bash
# 1. Configurar
cd pipeline_metagenomico
source setup_tmp_part4.sh

# 2. Copiar datos
cp /ruta/reads/*.fastq.gz data/MergedFastq/

# 3. Ejecutar
bash metagenomics_pipeline.sh
# Seleccionar: A (todos)
# Configurar Kraken2: 2 (Dual)
# Configurar bins: 2 (DAS Tool)
# Ejecutar: E
```

### Ejemplo 2: Solo Binning y Taxonomía

```bash
# Ejecutar módulos 3-6
bash metagenomics_pipeline.sh
# Seleccionar: 3, 4, 5, 6
# Ejecutar: E
```

### Ejemplo 3: Comparar Binnners

```bash
# Ejecutar GTDB-Tk con diferentes bins
export BINS_SOURCE=dastool
bash run_gtdbtk.sh

export BINS_SOURCE=metabat2
OUTPUT_DIR=output/gtdbtk_metabat2 bash run_gtdbtk.sh

export BINS_SOURCE=maxbin2
OUTPUT_DIR=output/gtdbtk_maxbin2 bash run_gtdbtk.sh
```

---

## 🐛 Solución de Problemas

### Error: "No space left on device"

**Causa:** Disco lleno, especialmente `/tmp` (root).

**Solución:**
```bash
# Configurar temporales en disco con espacio
source setup_tmp_part4.sh
echo $TMPDIR  # Verificar que apunta a disco con espacio
```

**[Ver guía completa →](ANALISIS_ARCHIVOS_TEMPORALES.md)**

### Error: GTDB-Tk "ValueError: __StageLogger"

**Causa:** Incompatibilidad con Python 3.14.

**Solución:**
```bash
# Recrear ambiente con Python 3.10
micromamba env remove -n gtdbtk
micromamba create -n gtdbtk -c bioconda python=3.10 gtdbtk=2.5.2
```

**[Ver solución completa →](SOLUCION_ERROR_GTDBTK.md)**

### Error: MaxBin2 "Failed to get Abundance information"

**Causa:** Formato incorrecto de archivo de abundancia.

**Solución:** El script actualizado usa reads directamente (`-reads`) en lugar de abundancia (`-abund`).

**[Ver corrección →](NOTA_CORRECCION_BINNING.md)**

---

## 🤝 Contribuciones

¡Las contribuciones son bienvenidas! Por favor:

1. Fork el repositorio
2. Crea una rama para tu feature (`git checkout -b feature/AmazingFeature`)
3. Commit tus cambios (`git commit -m 'Add some AmazingFeature'`)
4. Push a la rama (`git push origin feature/AmazingFeature`)
5. Abre un Pull Request

---

## 📄 Licencia

Este proyecto está bajo la Licencia MIT. Ver el archivo [LICENSE](LICENSE) para más detalles.

---

## 📧 Contacto

**Shaday Guerrero**

- GitHub: [@shadayguerrero](https://github.com/shadayguerrero)
- Email: shaday.guerrero@cinvestav.mx

**Repositorio:** [https://github.com/shadayguerrero/pipeline_metagenomico](https://github.com/shadayguerrero/pipeline_metagenomico)

---

## 🙏 Agradecimientos

Este pipeline integra las siguientes herramientas de código abierto:

- [Trim Galore](https://github.com/FelixKrueger/TrimGalore) - QC y trimming
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) - Mapeo de reads
- [MEGAHIT](https://github.com/voutcn/megahit) - Ensamblaje metagenómico
- [MetaBAT2](https://bitbucket.org/berkeleylab/metabat) - Binning
- [MaxBin2](https://sourceforge.net/projects/maxbin2/) - Binning
- [CONCOCT](https://github.com/BinPro/CONCOCT) - Binning
- [DAS Tool](https://github.com/cmks/DAS_Tool) - Refinamiento de bins
- [Kraken2](https://github.com/DerrickWood/kraken2) - Clasificación taxonómica
- [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) - Taxonomía de genomas
- [Prokka](https://github.com/tseemann/prokka) - Anotación de genes
- [RGI](https://github.com/arpcard/rgi) - Genes de resistencia
- [AntiSMASH](https://github.com/antismash/antismash) - Metabolitos secundarios

---

## 📊 Citación

Si usas este pipeline en tu investigación, por favor cita:

```bibtex
@software{guerrero2025pipeline,
  author = {Guerrero, Shaday},
  title = {Pipeline Metagenómico Modular},
  year = {2025},
  url = {https://github.com/shadayguerrero/pipeline_metagenomico}
}
```

Y las herramientas individuales según corresponda.

---

## 📈 Estadísticas

![GitHub stars](https://img.shields.io/github/stars/shadayguerrero/pipeline_metagenomico?style=social)
![GitHub forks](https://img.shields.io/github/forks/shadayguerrero/pipeline_metagenomico?style=social)
![GitHub issues](https://img.shields.io/github/issues/shadayguerrero/pipeline_metagenomico)

---

**Última actualización:** Diciembre 2025  
**Versión:** 1.0  
**Estado:** Activo y mantenido
