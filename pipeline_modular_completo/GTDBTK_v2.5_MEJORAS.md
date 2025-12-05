# GTDB-Tk v2.5.2 - Mejoras y Configuración

## Novedades en Esta Versión

El script de GTDB-Tk ha sido actualizado para aprovechar las nuevas características de la versión 2.5.2 con la base de datos r226.

---

## 🚀 Mejoras Principales

### 1. Soporte para Skani ✨ NUEVO

**¿Qué es Skani?**
- Herramienta ultrarrápida para calcular ANI (Average Nucleotide Identity)
- Reemplaza a FastANI en GTDB-Tk v2.5+
- **10-100x más rápido** que FastANI

**Configuración automática:**
```bash
# El script detecta automáticamente si tu base de datos tiene Skani
if [ -d "${GTDBTK_DATA_PATH}/skani" ]; then
    export GTDBTK_DISABLE_SKANI=0  # Habilitado
else
    export GTDBTK_DISABLE_SKANI=1  # Deshabilitado
fi
```

**Tu configuración:**
- Base de datos: `/data/database/gtdbtk_251103/gtdbtk_data_release226`
- Carpeta Skani: ✅ Presente (29,761 archivos)
- Estado: **Skani habilitado** (clasificación más rápida)

---

### 2. Selección Inteligente de Bins 🎯

El script ahora busca bins en **orden de preferencia**:

```
1. DAS Tool (refinados)     ← RECOMENDADO
   └─ dastool/DASTool_DASTool_bins/*.fa
   
2. MetaBAT2
   └─ metabat2/bin.*.fa
   
3. MaxBin2
   └─ maxbin2/bin.*.fasta
   
4. CONCOCT
   └─ concoct/bins/*.fa
```

**Ventaja:** Usa automáticamente los mejores bins disponibles sin configuración manual.

---

### 3. Base de Datos r226 (Actualizada)

**Cambios respecto a r207:**
- **+15,000 genomas** nuevos
- Taxonomía mejorada para grupos difíciles
- Mejor resolución a nivel de especie
- Soporte nativo para Skani

**Archeas:** Ahora usa `ar53` en lugar de `ar122`
- Antes: `gtdbtk.ar122.summary.tsv`
- Ahora: `gtdbtk.ar53.summary.tsv`

---

### 4. Mejor Manejo de Errores

**Verificaciones automáticas:**
- ✅ GTDB-Tk instalado y en PATH
- ✅ Base de datos existe
- ✅ Espacio en disco suficiente (>150 GB)
- ✅ Bins disponibles
- ✅ Versión de GTDB-Tk compatible

**Mensajes de error claros:**
```bash
✗ GTDB-Tk no está disponible

Solución:
  1. Activa el ambiente: micromamba activate gtdbtk
  2. Verifica la instalación: gtdbtk --version
```

---

### 5. Resumen Combinado Mejorado

**Archivos generados:**
```
output/gtdbtk/
├── GTDBTK_All_Bacteria.tsv    # Todas las bacterias
├── GTDBTK_All_Archaea.tsv     # Todas las arqueas
│
└── ${SAMPLE}/
    ├── gtdbtk.bac120.summary.tsv
    ├── gtdbtk.ar53.summary.tsv
    ├── gtdbtk.log
    └── ...
```

**Formato del resumen:**
- Encabezado preservado
- Todas las muestras combinadas
- Fácil de importar en R/Python

---

## 📊 Comparación de Versiones

| Característica | v2.4 (r207) | v2.5.2 (r226) |
|----------------|-------------|---------------|
| **Genomas en DB** | ~318,000 | ~333,000 |
| **Skani** | ❌ No | ✅ Sí |
| **Velocidad ANI** | Lento (FastANI) | Rápido (Skani) |
| **Arqueas** | ar122 | ar53 |
| **Taxonomía** | r207 | r226 (actualizada) |
| **Tiempo típico** | 4-8 horas | 2-4 horas |

---

## 🔧 Configuración

### Variables de Entorno

```bash
# Base de datos (REQUERIDO)
export GTDBTK_DATA_PATH="/data/database/gtdbtk_251103/gtdbtk_data_release226"

# Skani (OPCIONAL - detección automática)
export GTDBTK_DISABLE_SKANI=0  # 0=habilitado, 1=deshabilitado

# Directorio temporal (OPCIONAL)
export TMPDIR="/tmp/gtdbtk"
```

### Parámetros del Script

```bash
# Directorios
BINNING_DIR="${BINNING_DIR:-output/binning}"
OUTPUT_DIR="${OUTPUT_DIR:-output/gtdbtk}"

# Recursos
THREADS="${THREADS:-40}"
```

---

## 🚀 Uso

### Opción 1: Ejecutar Directamente

```bash
# Activar ambiente
micromamba activate gtdbtk

# Configurar base de datos
export GTDBTK_DATA_PATH="/data/database/gtdbtk_251103/gtdbtk_data_release226"

# Ejecutar
bash run_gtdbtk.sh
```

### Opción 2: Desde el Pipeline Principal

```bash
# Ejecutar pipeline
bash metagenomics_pipeline.sh

# Seleccionar módulo 6 (GTDB-Tk)
```

### Opción 3: Con Variables Personalizadas

```bash
# Activar ambiente
micromamba activate gtdbtk

# Ejecutar con configuración personalizada
BINNING_DIR="/ruta/custom/binning" \
OUTPUT_DIR="/ruta/custom/gtdbtk" \
THREADS=60 \
bash run_gtdbtk.sh
```

---

## 📈 Salida Esperada

### Ejemplo de Ejecución

```
========================================
GTDB-Tk - Clasificación Taxonómica
========================================

Configuración:
  Base de datos: /data/database/gtdbtk_251103/gtdbtk_data_release226
  Directorio de bins: output/binning
  Directorio de salida: output/gtdbtk
  Threads: 40
  Skani: ✓ Habilitado (rápido)
  TMPDIR: /tmp/gtdbtk_12345

    ✓ Espacio temporal disponible: 250 GB
    ✓ GTDB-Tk v2.5.2 detectado
    ✓ Base de datos: release226

Muestras encontradas:
  - SRR5936076
  - SRR5936077

========================================
Procesando: SRR5936076
========================================
  Fuente de bins: DAS Tool (refinados)
  Bins encontrados: 12
  Directorio: output/binning/SRR5936076/dastool/DASTool_DASTool_bins

  Ejecutando GTDB-Tk classify_wf...

    ✓ Clasificación completada

  Resultados:
    Bacterias: 11 bins
    Arqueas: 1 bins

========================================
Generando Resumen Combinado
========================================

  Combinando clasificaciones bacterianas...
    ✓ 22 bins bacterianos clasificados
  Combinando clasificaciones arqueales...
    ✓ 2 bins arqueales clasificados

========================================
GTDB-Tk COMPLETADO
========================================

Resultados:
  ✓ Exitosas: 2
  ✗ Fallidas: 0

Clasificaciones:
  Bacterias: 22 bins → output/gtdbtk/GTDBTK_All_Bacteria.tsv
  Arqueas: 2 bins → output/gtdbtk/GTDBTK_All_Archaea.tsv

Archivos individuales en: output/gtdbtk/
```

---

## 🔍 Interpretación de Resultados

### Archivo de Resumen (gtdbtk.bac120.summary.tsv)

**Columnas principales:**
- `user_genome` - Nombre del bin
- `classification` - Taxonomía completa (d__; p__; c__; o__; f__; g__; s__)
- `fastani_reference` - Genoma de referencia más cercano
- `fastani_ani` - ANI con la referencia (%)
- `classification_method` - Método usado (taxonomy, ANI, etc.)

**Ejemplo:**
```
user_genome    classification                                      fastani_ani
bin_1.fa       d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__Pseudomonas aeruginosa    98.5
```

### Niveles Taxonómicos

- `d__` - Dominio (Domain)
- `p__` - Filo (Phylum)
- `c__` - Clase (Class)
- `o__` - Orden (Order)
- `f__` - Familia (Family)
- `g__` - Género (Genus)
- `s__` - Especie (Species)

---

## ⏱️ Tiempo de Ejecución

**Con Skani (r226):**
- 10 bins: 20-40 min
- 50 bins: 1-2 horas
- 100 bins: 2-4 horas

**Sin Skani (r207):**
- 10 bins: 1-2 horas
- 50 bins: 4-8 horas
- 100 bins: 8-16 horas

**Mejora:** ~3-4x más rápido con Skani

---

## 🐛 Solución de Problemas

### Error: "GTDBTK_DATA_PATH not set"

```bash
export GTDBTK_DATA_PATH="/data/database/gtdbtk_251103/gtdbtk_data_release226"
```

### Error: "No bins found"

Verifica que el binning se completó:
```bash
ls output/binning/*/dastool/DASTool_DASTool_bins/*.fa
```

### Error: "Out of memory"

Reduce threads o aumenta RAM:
```bash
THREADS=20 bash run_gtdbtk.sh
```

### Advertencia: "Skani disabled"

Tu base de datos no tiene la carpeta `skani/`. Opciones:
1. Actualizar a r226 completa
2. Continuar sin Skani (más lento pero funcional)

---

## 📚 Referencias

- **GTDB-Tk:** Chaumeil et al. (2022) Bioinformatics 38(23):5315-5316
- **GTDB r226:** Parks et al. (2022) Nucleic Acids Research 50(D1):D785-D794
- **Skani:** Jain et al. (2023) bioRxiv

---

**Última actualización:** 2 de diciembre de 2025  
**Versión del script:** 2.0 (compatible con GTDB-Tk v2.5.2 y r226)
