# MetaConexus: Modular Metagenomic Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Platform](https://img.shields.io/badge/Platform-Linux-blue.svg)](https://www.linux.org/)
[![Shell](https://img.shields.io/badge/Shell-Bash-green.svg)](https://www.gnu.org/software/bash/)
[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/)

Complete and modular pipeline for high-throughput metagenomic analysis, from raw reads to integrated final reports.

---

## 📋 Table of Contents

- [Description](#description)
- [Features](#features)
- [Workflow Diagram](#workflow-diagram)
- [Pipeline Modules](#pipeline-modules)
- [Quick Installation](#quick-installation)
- [Usage](#usage)
- [Documentation](#documentation)
- [Requirements](#requirements)
- [Directory Structure](#directory-structure)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)
- [Contributions](#contributions)
- [License](#license)
- [Contact](#contact)

---

## 🧬 Description

This metagenomic pipeline processes next-generation sequencing (NGS) data to characterize complex microbial communities. It integrates state-of-the-art bioinformatics tools to:

- **Preprocess** sequencing reads
- **Assemble** metagenomic genomes
- **Classify** reads and genomes taxonomically
- **Annotate** genes and proteins functionally
- **Identify** antibiotic resistance genes
- **Detect** biosynthetic gene clusters for secondary metabolites
- **Generate** integrated reports and visualizations

Designed to be **modular**, **scalable**, and **user-friendly**, allowing execution of the complete pipeline or individual modules according to project needs.

---

## ✨ Features

- ✅ **Modular:** Run individual modules or the complete pipeline
- ✅ **Flexible:** Bin selection (DAS Tool, MetaBAT2, MaxBin2, CONCOCT)
- ✅ **Robust:** Binning with 4 tools + refinement with DAS Tool
- ✅ **Updated:** GTDB-Tk v2.5.2 with database r226 and Skani
- ✅ **Complete:** 10 modules from QC to integrative analysis
- ✅ **Optimized:** Intelligent temporary file management
- ✅ **Documented:** Detailed guides for each step
- ✅ **Reproducible:** Versioned scripts and explicit configuration

---

## 📊 Workflow Diagram

![Pipeline Workflow Diagram](https://private-us-east-1.manuscdn.com/sessionFile/qaU4ooowr2qViGx42i6Pxe/sandbox/05XdRNi9wICX0Gk2LbQtpz-images_1765309228082_na1fn_L2hvbWUvdWJ1bnR1L3BpcGVsaW5lX2Zsb3c.png?Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9wcml2YXRlLXVzLWVhc3QtMS5tYW51c2Nkbi5jb20vc2Vzc2lvbkZpbGUvcWFVNG9vb3dyMnFWaUd4NDJpNlB4ZS9zYW5kYm94LzA1WGRSTmk5d0lDWDBHazJMYlF0cHotaW1hZ2VzXzE3NjUzMDkyMjgwODJfbmExZm5fTDJodmJXVXZkV0oxYm5SMUwzQnBjR1ZzYVc1bFgyWnNiM2MucG5nIiwiQ29uZGl0aW9uIjp7IkRhdGVMZXNzVGhhbiI6eyJBV1M6RXBvY2hUaW1lIjoxNzk4NzYxNjAwfX19XX0_&Key-Pair-Id=K2HSFNDJXOU9YS&Signature=eEps~5Tl~hKPhRIguS5U1TZpXyXX4FChbdAcTCG~tQhgA8oe5JLmUZebPso0nzWJbjuwGcEBIAJvY09k~xnredgoD-H56Zya5V3L0~GWXSPnD5rHRdgFM1lY2j1i-OU5XGtn8cK8lNmJ5aLY~Ibf3toAS-lTRZpNlYKvIBIWnIeX1WXFDTyqX0laICIn54qNF0mQmvg7J1qntSF-Xt59PQ0iQPTzaUzZZC7y0LclY3WH4OVMdiGabEOFLLeTUaWy2mjAEN~Kyc0rnMYpxnQLOit7yQziS7oTmMRmjhK-aVvYB2fvtMFwMYto9qUwJtobz5rqDr7wPwryi~GnjgsKvA__)

**[View complete workflow documentation →](https://github.com/shadayguerrero/pipeline_metagenomico/blob/main/pipeline_modular_completo/FLUJO_PIPELINE.md)**

---

## 🔧 Pipeline Modules

| # | Module | Tool | Function | Time* |
|---|--------|------|----------|-------|
| **1** | QC & Trimming | Trim Galore | Quality control and read cleaning | 30-60 min |
| **2** | Host Removal | Bowtie2 | Host DNA removal | 1-2 h |
| **3** | Assembly | MEGAHIT | De novo contig assembly | 4-8 h |
| **4** | Binning | MetaBAT2, MaxBin2, CONCOCT, DAS Tool | Genome reconstruction (MAGs) | 2-3 h |
| **5** | Read Taxonomy | Kraken2 | Taxonomic classification of reads | 30-60 min |
| **6** | Bin Taxonomy | GTDB-Tk v2.5.2 | Taxonomic classification of genomes | 2-4 h |
| **7** | Annotation | Prokka | Functional gene annotation | 1-2 h |
| **8** | Resistance | RGI + CARD | Antibiotic resistance gene identification | 30-60 min |
| **9** | Metabolites | AntiSMASH | Biosynthetic cluster detection | 2-4 h |
| **10** | Analysis | Python + R | Integrated reports and visualizations | 30-60 min |

\* *Estimated time for 2 samples with 40 threads*

**Total time:** 23-37 hours (~1-1.5 days)

---

## 🚀 Quick Installation

### Prerequisites

- Linux (Ubuntu 20.04+ or CentOS 7+)
- [Micromamba](https://mamba.readthedocs.io/en/latest/installation.html) or Conda
- 40+ threads recommended
- 500+ GB disk space

### Step 1: Clone the Repository

```bash
git clone https://github.com/shadayguerrero/pipeline_metagenomico.git
cd pipeline_metagenomico
```

### Step 2: Install Environments

```bash
# Create environments with micromamba
bash INSTALACION_AMBIENTES.sh
```

### Step 3: Configure Databases

```bash
# GTDB-Tk (required for module 6)
export GTDBTK_DATA_PATH="/path/to/gtdbtk_data_release226"

# Kraken2 (required for module 5)
# Download databases from https://benlangmead.github.io/aws-indexes/k2
```

### Step 4: Configure Temporary Directories

```bash
# If your root is full, configure temporaries on disk with space
source setup_tmp_part4.sh
```

**[View complete configuration guide →](https://github.com/shadayguerrero/pipeline_metagenomico/blob/main/pipeline_modular_completo/CONFIGURACION_RAPIDA.md)**

---

## 💻 Usage

### Interactive Execution (Recommended)

```bash
# Run the pipeline with interactive menu
bash metagenomics_pipeline.sh
```

**Main menu:**
- Select individual modules (1-10) or all (A)
- Configure Kraken2 mode (Simple, Dual, Triple)
- Select bin source (DAS Tool, MetaBAT2, etc.)
- Review selection (R) and execute (E)

### Individual Module Execution

```bash
# Module 1: QC & Trimming
bash run_trimgalore.sh

# Module 4: Binning
bash run_binning_fixed.sh

# Module 6: GTDB-Tk
export BINS_SOURCE=dastool
bash run_gtdbtk.sh
```

### Complete Automated Execution

```bash
# Configure variables
export BINS_SOURCE=dastool
export KRAKEN_MODE=dual

# Run complete pipeline
bash metagenomics_pipeline.sh
# Select: A (all modules)
# Execute: E
```

**[View detailed usage guide →](GUIA_RAPIDA.md)**

---

## 📚 Documentation

### Getting Started Guides

- **[Quick Configuration](CONFIGURACION_RAPIDA.md)** - Installation and configuration in 5 minutes
- **[Usage Guide](GUIA_RAPIDA.md)** - Use cases and practical examples
- **[Pipeline Workflow](FLUJO_PIPELINE.md)** - Detailed description of each module

### Technical Documentation

- **[Bin Selection](SELECCION_BINS.md)** - How to choose between DAS Tool, MetaBAT2, MaxBin2, CONCOCT
- **[Complete Binning](BINNING_COMPLETO.md)** - Binning with 4 tools and refinement
- **[GTDB-Tk v2.5.2](GTDBTK_v2.5_MEJORAS.md)** - GTDB-Tk updates and database r226
- **[Temporary Files](ANALISIS_ARCHIVOS_TEMPORALES.md)** - Disk space management
- **[Usage on /Part4](GUIA_USO_PART4.md)** - Configuration for external disks

### Troubleshooting

- **[Binning Correction](NOTA_CORRECCION_BINNING.md)** - Common binning errors
- **[GTDB-Tk Error](SOLUCION_ERROR_GTDBTK.md)** - Python 3.14 incompatibility

---

## 🛠️ Requirements

### Hardware

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| **CPU** | 20 cores | 40+ cores |
| **RAM** | 64 GB | 128+ GB |
| **Disk** | 500 GB | 1+ TB |
| **Temporary** | 300 GB | 500+ GB |

### Software

- **Operating System:** Linux (Ubuntu 20.04+, CentOS 7+)
- **Environment manager:** Micromamba or Conda
- **Python:** 3.8-3.11 (NOT 3.12+)
- **Bash:** 4.0+

### Databases

| Database | Size | Module | Download |
|----------|------|--------|----------|
| **GTDB-Tk r226** | ~80 GB | 6 | [Link](https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_r226_data.tar.gz) |
| **Kraken2 GTDB** | ~60 GB | 5 | [Link](https://benlangmead.github.io/aws-indexes/k2) |
| **CARD** | ~1 GB | 8 | [Link](https://card.mcmaster.ca/download) |
| **AntiSMASH** | ~15 GB | 9 | Installed with environment |

---

## 📁 Directory Structure

```
pipeline_metagenomico/
├── metagenomics_pipeline.sh       # Main pipeline
├── run_*.sh                       # Module scripts
├── *.yaml                         # Environment files
├── setup_tmp_part4.sh             # Temporary configuration
│
├── docs/                          # Documentation
│   ├── CONFIGURACION_RAPIDA.md
│   ├── FLUJO_PIPELINE.md
│   ├── GUIA_RAPIDA.md
│   └── ...
│
└── data/                          # Data (create manually)
    ├── MergedFastq/               # Input reads
    ├── output/                    # Pipeline results
    └── tmp/                       # Temporary files
```

### Pipeline Outputs

```
output/
├── trimmed/                       # Module 1: Clean reads
├── host_removed/                  # Module 2: Host-free reads
├── megahit_assemblies/            # Module 3: Assembled contigs
├── binning/                       # Module 4: Reconstructed genomes
│   ├── metabat2/
│   ├── maxbin2/
│   ├── concoct/
│   └── dastool/                   # ← Refined bins (recommended)
├── kraken2/                       # Module 5: Read taxonomy
├── gtdbtk/                        # Module 6: Bin taxonomy
├── prokka/                        # Module 7: Annotated genes
├── rgi/                           # Module 8: Resistance genes
├── antismash/                     # Module 9: Biosynthetic clusters
└── analysis/                      # Module 10: Final reports
```

---

## 💡 Examples

### Example 1: Complete Pipeline

```bash
# 1. Configure
cd pipeline_metagenomico
source setup_tmp_part4.sh

# 2. Copy data
cp /path/reads/*.fastq.gz data/MergedFastq/

# 3. Execute
bash metagenomics_pipeline.sh
# Select: A (all)
# Configure Kraken2: 2 (Dual)
# Configure bins: 2 (DAS Tool)
# Execute: E
```

### Example 2: Binning and Taxonomy Only

```bash
# Run modules 3-6
bash metagenomics_pipeline.sh
# Select: 3, 4, 5, 6
# Execute: E
```

### Example 3: Compare Binners

```bash
# Run GTDB-Tk with different bins
export BINS_SOURCE=dastool
bash run_gtdbtk.sh

export BINS_SOURCE=metabat2
OUTPUT_DIR=output/gtdbtk_metabat2 bash run_gtdbtk.sh

export BINS_SOURCE=maxbin2
OUTPUT_DIR=output/gtdbtk_maxbin2 bash run_gtdbtk.sh
```

---

## 🐛 Troubleshooting

### Error: "No space left on device"

**Cause:** Disk full, especially `/tmp` (root).

**Solution:**
```bash
# Configure temporaries on disk with space
source setup_tmp_part4.sh
echo $TMPDIR  # Verify it points to disk with space
```

**[View complete guide →](ANALISIS_ARCHIVOS_TEMPORALES.md)**

### Error: GTDB-Tk "ValueError: __StageLogger"

**Cause:** Incompatibility with Python 3.14.

**Solution:**
```bash
# Recreate environment with Python 3.10
micromamba env remove -n gtdbtk
micromamba create -n gtdbtk -c bioconda python=3.10 gtdbtk=2.5.2
```

**[View complete solution →](SOLUCION_ERROR_GTDBTK.md)**

### Error: MaxBin2 "Failed to get Abundance information"

**Cause:** Incorrect abundance file format.

**Solution:** The updated script uses reads directly (`-reads`) instead of abundance (`-abund`).

**[View correction →](NOTA_CORRECCION_BINNING.md)**

---

## 🤝 Contributions

Contributions are welcome! Please:

1. Fork the repository
2. Create a branch for your feature (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

---

## 📄 Acknowledgments

This project was carried out during my Visiting Professor stay at CINVESTAV Unidad Irapuato, in the Bioinformatics and Complex Networks Laboratory directed by Dr. Maribel Hernandez. The work was done in collaboration with Dr. Elizabeth Cadenas (elizabeth.cadenas.c@gmail.com), Dr. Katia Aviña-Padilla (Katia.avinap@cinvestav.mx), for the study of Las Pozas de Cuatro Ciénegas, Coahuila.

---

## 📧 Contact

**Shaday Guerrero**

- GitHub: [@shadayguerrero](https://github.com/shadayguerrero)
- Email: shaday@matmor.unam.mx

**Repository:** [https://github.com/shadayguerrero/pipeline_metagenomico](https://github.com/shadayguerrero/pipeline_metagenomico)

---

## 🙏 Credits

This pipeline integrates the following open-source tools:

- [Trim Galore](https://github.com/FelixKrueger/TrimGalore) - QC and trimming
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) - Read mapping
- [MEGAHIT](https://github.com/voutcn/megahit) - Metagenomic assembly
- [MetaBAT2](https://bitbucket.org/berkeleylab/metabat) - Binning
- [MaxBin2](https://sourceforge.net/projects/maxbin2/) - Binning
- [CONCOCT](https://github.com/BinPro/CONCOCT) - Binning
- [DAS Tool](https://github.com/cmks/DAS_Tool) - Bin refinement
- [Kraken2](https://github.com/DerrickWood/kraken2) - Taxonomic classification
- [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) - Genome taxonomy
- [Prokka](https://github.com/tseemann/prokka) - Gene annotation
- [RGI](https://github.com/arpcard/rgi) - Resistance genes
- [AntiSMASH](https://github.com/antismash/antismash) - Secondary metabolites

---

And the individual tools as appropriate.

---

## 📈 Statistics

![GitHub stars](https://img.shields.io/github/stars/shadayguerrero/pipeline_metagenomico?style=social)
![GitHub forks](https://img.shields.io/github/forks/shadayguerrero/pipeline_metagenomico?style=social)
![GitHub issues](https://img.shields.io/github/issues/shadayguerrero/pipeline_metagenomico)

---

**Last updated:** December 2025  
**Version:** 1.0  
**Status:** Active and maintained
