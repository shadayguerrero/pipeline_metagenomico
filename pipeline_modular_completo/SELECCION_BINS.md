# Selección de Bins para Módulos 6-9

Los módulos de análisis posterior (GTDB-Tk, Prokka, RGI y AntiSMASH) ahora permiten seleccionar qué conjunto de bins usar mediante la variable de entorno `BINS_SOURCE`.

---

## 🎯 Opciones Disponibles

### 1. **auto** (Por defecto)
Busca bins en orden de preferencia:
1. DAS Tool (refinados) ← Recomendado
2. MetaBAT2
3. MaxBin2
4. CONCOCT

**Uso:**
```bash
bash run_gtdbtk.sh
# O explícitamente:
BINS_SOURCE=auto bash run_gtdbtk.sh
```

---

### 2. **dastool**
Usa solo bins refinados de DAS Tool.

**Uso:**
```bash
BINS_SOURCE=dastool bash run_gtdbtk.sh
BINS_SOURCE=dastool bash run_prokka.sh
BINS_SOURCE=dastool bash run_rgi.sh
BINS_SOURCE=dastool bash run_antismash.sh
```

**Recomendado para:**
- Análisis final
- Máxima calidad de bins
- Cuando quieres los mejores bins de los 3 binnners

---

### 3. **metabat2**
Usa solo bins de MetaBAT2.

**Uso:**
```bash
BINS_SOURCE=metabat2 bash run_gtdbtk.sh
```

**Recomendado para:**
- Comparar con otros binnners
- Cuando DAS Tool no generó bins
- Bins de alta cobertura

---

### 4. **maxbin2**
Usa solo bins de MaxBin2.

**Uso:**
```bash
BINS_SOURCE=maxbin2 bash run_gtdbtk.sh
```

**Recomendado para:**
- Comparar con otros binnners
- Bins con buenos marcadores de copia única
- Bacterias y arqueas bien caracterizadas

---

### 5. **concoct**
Usa solo bins de CONCOCT.

**Uso:**
```bash
BINS_SOURCE=concoct bash run_gtdbtk.sh
```

**Recomendado para:**
- Comparar con otros binnners
- Genomas de baja abundancia
- Cuando quieres más bins (CONCOCT genera muchos)

---

## 📊 Comparación de Binnners

| Característica | DAS Tool | MetaBAT2 | MaxBin2 | CONCOCT |
|----------------|----------|----------|---------|---------|
| **Tipo** | Refinamiento | Binning | Binning | Binning |
| **Método** | Combina los 3 | Cobertura + composición | Marcadores + cobertura | Fragmentación + clustering |
| **Bins típicos** | 5-30 | 10-50 | 10-40 | 20-100 |
| **Calidad** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐ |
| **Completitud** | Alta | Alta | Media-Alta | Media |
| **Contaminación** | Baja | Baja-Media | Baja | Media |
| **Mejor para** | Análisis final | Genomas abundantes | Bacterias/arqueas | Genomas raros |

---

## 🔄 Ejecutar con Diferentes Bins

### Ejemplo 1: Comparar GTDB-Tk con 3 Binnners

```bash
# Activar ambiente
micromamba activate gtdbtk

# DAS Tool (refinados)
BINS_SOURCE=dastool OUTPUT_DIR=output/gtdbtk_dastool bash run_gtdbtk.sh

# MetaBAT2
BINS_SOURCE=metabat2 OUTPUT_DIR=output/gtdbtk_metabat2 bash run_gtdbtk.sh

# MaxBin2
BINS_SOURCE=maxbin2 OUTPUT_DIR=output/gtdbtk_maxbin2 bash run_gtdbtk.sh
```

### Ejemplo 2: Pipeline Completo con DAS Tool

```bash
# Módulo 6: GTDB-Tk
micromamba activate gtdbtk
BINS_SOURCE=dastool bash run_gtdbtk.sh

# Módulo 7: Prokka
micromamba activate prokka
BINS_SOURCE=dastool bash run_prokka.sh

# Módulo 8: RGI
micromamba activate rgi
BINS_SOURCE=dastool bash run_rgi.sh

# Módulo 9: AntiSMASH
micromamba activate antismash
BINS_SOURCE=dastool bash run_antismash.sh
```

### Ejemplo 3: Pipeline Completo con MetaBAT2

```bash
# Usar MetaBAT2 para todos los módulos
export BINS_SOURCE=metabat2

micromamba activate gtdbtk && bash run_gtdbtk.sh
micromamba activate prokka && bash run_prokka.sh
micromamba activate rgi && bash run_rgi.sh
micromamba activate antismash && bash run_antismash.sh
```

---

## 📁 Estructura de Directorios

```
output/binning/SAMPLE/
├── metabat2/
│   ├── bin.1.fa
│   ├── bin.2.fa
│   └── ...
│
├── maxbin2/
│   ├── bin.001.fasta
│   ├── bin.002.fasta
│   └── ...
│
├── concoct/
│   └── bins/
│       ├── 1.fa
│       ├── 2.fa
│       └── ...
│
└── dastool/
    └── DASTool_DASTool_bins/  ← RECOMENDADO
        ├── bin_1.fa
        ├── bin_2.fa
        └── ...
```

---

## ⚙️ Variables de Entorno

Todas las variables disponibles para personalizar:

```bash
# Selección de bins
export BINS_SOURCE=dastool  # auto, dastool, metabat2, maxbin2, concoct

# Directorios
export BINNING_DIR=/ruta/custom/binning
export OUTPUT_DIR=/ruta/custom/salida

# Recursos
export THREADS=40

# GTDB-Tk específico
export GTDBTK_DATA_PATH=/ruta/a/gtdbtk_db
export TMPDIR=/tmp/gtdbtk

# Ejecutar
bash run_gtdbtk.sh
```

---

## 🔍 Verificar Qué Bins Se Usaron

Cada script muestra en la salida qué fuente de bins está usando:

```
========================================
GTDB-Tk - Clasificación Taxonómica
========================================

Configuración:
  Base de datos: /data/database/gtdbtk_251103/gtdbtk_data_release226
  Directorio de bins: output/binning
  Directorio de salida: output/gtdbtk
  Fuente de bins: dastool         ← AQUÍ
  Threads: 40
  Skani: ✓ Habilitado (rápido)
  TMPDIR: /tmp/gtdbtk_12345

========================================
Procesando: SRR5936076
========================================
  Fuente de bins: DAS Tool (refinados)  ← Y AQUÍ
  Bins encontrados: 12
  Directorio: output/binning/SRR5936076/dastool/DASTool_DASTool_bins
```

---

## 💡 Recomendaciones

### Para Análisis Exploratorio
```bash
# Usa MetaBAT2 (rápido, buenos resultados)
BINS_SOURCE=metabat2 bash run_gtdbtk.sh
```

### Para Análisis Final
```bash
# Usa DAS Tool (máxima calidad)
BINS_SOURCE=dastool bash run_gtdbtk.sh
BINS_SOURCE=dastool bash run_prokka.sh
BINS_SOURCE=dastool bash run_rgi.sh
BINS_SOURCE=dastool bash run_antismash.sh
```

### Para Comparación
```bash
# Ejecuta con los 3 binnners y compara resultados
for binner in metabat2 maxbin2 concoct; do
    BINS_SOURCE=$binner \
    OUTPUT_DIR=output/gtdbtk_$binner \
    bash run_gtdbtk.sh
done
```

---

## ❓ Solución de Problemas

### Error: "No se encontraron bins de dastool"

**Causa:** DAS Tool no generó bins o falló.

**Solución:**
```bash
# Verificar si DAS Tool generó bins
ls output/binning/*/dastool/DASTool_DASTool_bins/*.fa

# Si no hay bins, usa otro binner
BINS_SOURCE=metabat2 bash run_gtdbtk.sh
```

### Error: "No se encontraron bins para SAMPLE"

**Causa:** El binner seleccionado no tiene bins para esa muestra.

**Solución:**
```bash
# Verificar qué binnners tienen bins
ls output/binning/SAMPLE/*/

# Usar modo auto para buscar automáticamente
BINS_SOURCE=auto bash run_gtdbtk.sh
```

---

## 📚 Referencias

- **MetaBAT2:** Kang et al. (2019) PeerJ 7:e7359
- **MaxBin2:** Wu et al. (2016) Bioinformatics 32(4):605-607
- **CONCOCT:** Alneberg et al. (2014) Nature Methods 11:1144-1146
- **DAS Tool:** Sieber et al. (2018) Nature Microbiology 3:836-843

---

**Última actualización:** 2 de diciembre de 2025  
**Versión:** 1.0
