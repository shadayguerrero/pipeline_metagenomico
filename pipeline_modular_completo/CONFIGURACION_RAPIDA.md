# Configuración Rápida para Nuevo Servidor

## 🎯 Objetivo

Configurar el pipeline metagenómico en el nuevo servidor usando `/mnt/Part4/Laboratory` para datos, resultados y archivos temporales.

---

## ⚡ Pasos Rápidos (5 minutos)

### 1. Transferir el Paquete

```bash
# En tu máquina local
scp pipeline_modular_completo.tar.gz usuario@servidor:/mnt/Part4/Laboratory/
```

### 2. Extraer en el Servidor

```bash
# En el servidor
cd /mnt/Part4/Laboratory
tar -xzf pipeline_modular_completo.tar.gz
cd pipeline_modular_completo
```

### 3. Crear Estructura de Directorios

```bash
# Crear directorios necesarios
mkdir -p /mnt/Part4/Laboratory/data/MergedFastq
mkdir -p /mnt/Part4/Laboratory/output
mkdir -p /mnt/Part4/Laboratory/tmp
```

### 4. Configurar Temporales

```bash
# Configurar directorios temporales
source setup_tmp_part4.sh
```

**Salida esperada:**
```
============================================================
  Configuración de Directorios Temporales
============================================================

Creando directorios temporales en /mnt/Part4/Laboratory/tmp...
  ✓ /mnt/Part4/Laboratory/tmp/general
  ✓ /mnt/Part4/Laboratory/tmp/megahit
  ✓ /mnt/Part4/Laboratory/tmp/gtdbtk
  ✓ /mnt/Part4/Laboratory/tmp/bowtie2
  ✓ /mnt/Part4/Laboratory/tmp/concoct
  ✓ /mnt/Part4/Laboratory/tmp/antismash
  ✓ /mnt/Part4/Laboratory/tmp/maxbin2
  ✓ /mnt/Part4/Laboratory/tmp/prokka

Configuración:
  Directorio base: /mnt/Part4/Laboratory/tmp
  Espacio disponible: 3.9T

Variables de entorno configuradas:
  TMPDIR=/mnt/Part4/Laboratory/tmp/general
  MEGAHIT_TMP=/mnt/Part4/Laboratory/tmp/megahit
  GTDBTK_TMP=/mnt/Part4/Laboratory/tmp/gtdbtk
  CONCOCT_TMP=/mnt/Part4/Laboratory/tmp/concoct
  ANTISMASH_TMP=/mnt/Part4/Laboratory/tmp/antismash

✓ Configuración completada
```

### 5. Editar Rutas en el Pipeline

```bash
# Editar el archivo principal
nano metagenomics_pipeline.sh
```

**Cambiar estas líneas (aproximadamente líneas 27-32):**

```bash
# ANTES (valores por defecto):
PROJECT_DIR="/files/shaday/4_cienegas"
OUTPUT_DIR="${PROJECT_DIR}/output"

# DESPUÉS (tu configuración):
PROJECT_DIR="/mnt/Part4/Laboratory/data"
OUTPUT_DIR="/mnt/Part4/Laboratory/output"
```

**Guardar:** `Ctrl+O`, `Enter`, `Ctrl+X`

### 6. Copiar Datos

```bash
# Copiar tus archivos FASTQ
cp /ruta/original/*.fastq.gz /mnt/Part4/Laboratory/data/MergedFastq/

# Verificar
ls -lh /mnt/Part4/Laboratory/data/MergedFastq/
```

### 7. Verificar Configuración

```bash
# Verificar espacio
df -h /mnt/Part4

# Verificar TMPDIR
echo $TMPDIR
# Debe mostrar: /mnt/Part4/Laboratory/tmp/general

# Verificar estructura
tree -L 2 /mnt/Part4/Laboratory
```

### 8. Ejecutar Pipeline

```bash
# Ejecutar
bash metagenomics_pipeline.sh
```

---

## 📁 Estructura Final

```
/mnt/Part4/Laboratory/
├── pipeline_modular_completo/
│   ├── metagenomics_pipeline.sh
│   ├── run_*.sh
│   ├── setup_tmp_part4.sh
│   └── *.yaml
│
├── data/
│   └── MergedFastq/
│       ├── SRR5936076_R1.fastq.gz
│       ├── SRR5936076_R2.fastq.gz
│       ├── SRR5936077_R1.fastq.gz
│       └── SRR5936077_R2.fastq.gz
│
├── output/
│   ├── trimmed/
│   ├── host_removed/
│   ├── megahit_assemblies/
│   ├── binning/
│   ├── gtdbtk/
│   ├── prokka/
│   ├── rgi/
│   └── antismash/
│
└── tmp/                    # SE PUEDE BORRAR DESPUÉS
    ├── general/
    ├── megahit/
    ├── gtdbtk/
    ├── concoct/
    └── antismash/
```

---

## 🔧 Variables de Entorno Importantes

Agregar al `~/.bashrc` o `~/.bash_profile`:

```bash
# Directorios temporales
export TMPDIR="/mnt/Part4/Laboratory/tmp/general"
export TEMP="${TMPDIR}"
export TMP="${TMPDIR}"

# GTDB-Tk (ajustar a tu ruta)
export GTDBTK_DATA_PATH="/data/database/gtdbtk_251103/gtdbtk_data_release226"

# Micromamba (si no está configurado)
export MAMBA_ROOT_PREFIX="/home/shaday/micromamba"
eval "$(micromamba shell hook --shell bash)"
```

Luego:

```bash
source ~/.bashrc
```

---

## ✅ Verificación Pre-Ejecución

```bash
# 1. Verificar espacio (>500 GB recomendado)
df -h /mnt/Part4
# Debe mostrar ~3.9 TiB disponible

# 2. Verificar TMPDIR
echo $TMPDIR
# Debe mostrar: /mnt/Part4/Laboratory/tmp/general

# 3. Verificar datos
ls /mnt/Part4/Laboratory/data/MergedFastq/*.fastq.gz
# Debe listar tus archivos

# 4. Verificar ambientes
micromamba env list
# Debe mostrar: qc_assembly, binning, kraken2, gtdbtk, prokka, rgi, antismash

# 5. Verificar GTDB-Tk
micromamba activate gtdbtk
gtdbtk check_install
# Debe mostrar: OK en todas las verificaciones

# 6. Verificar permisos
touch /mnt/Part4/Laboratory/tmp/test.txt && rm /mnt/Part4/Laboratory/tmp/test.txt
echo "✓ Permisos OK"
```

---

## 🚀 Ejecución

### Opción 1: Pipeline Completo

```bash
cd /mnt/Part4/Laboratory/pipeline_modular_completo
source setup_tmp_part4.sh
bash metagenomics_pipeline.sh
# Seleccionar: A (todos los módulos)
# Configurar Kraken2 si es necesario
# Configurar bins: 2 (DAS Tool)
# Ejecutar: E
```

### Opción 2: Módulos Individuales

```bash
# Ejecutar solo módulos específicos
bash metagenomics_pipeline.sh
# Seleccionar módulos individuales (1-10)
# Ejecutar: E
```

### Opción 3: Scripts Individuales

```bash
# Activar ambiente y ejecutar script
micromamba activate gtdbtk
export BINS_SOURCE=dastool
bash run_gtdbtk.sh
```

---

## 🧹 Limpieza Post-Ejecución

```bash
# Limpiar archivos temporales (libera ~320 GB)
rm -rf /mnt/Part4/Laboratory/tmp/*

# Recrear estructura
mkdir -p /mnt/Part4/Laboratory/tmp/{general,megahit,gtdbtk,concoct,antismash,maxbin2,prokka}

# NO borrar output/ (contiene resultados)
```

---

## 📊 Monitoreo Durante Ejecución

### Terminal 1: Pipeline

```bash
bash metagenomics_pipeline.sh
```

### Terminal 2: Monitoreo

```bash
# Ver uso de espacio
watch -n 10 'df -h /mnt/Part4'

# Ver archivos temporales
watch -n 30 'du -sh /mnt/Part4/Laboratory/tmp/*'

# Ver procesos
htop
```

---

## ⚠️ Solución de Problemas Comunes

### Error: "No space left on device"

```bash
# Verificar que TMPDIR está configurado
echo $TMPDIR
# Si está vacío:
source setup_tmp_part4.sh
```

### Error: "Permission denied"

```bash
# Dar permisos
sudo chmod -R 1777 /mnt/Part4/Laboratory/tmp
# O cambiar propietario
sudo chown -R $USER:$USER /mnt/Part4/Laboratory
```

### Error: "GTDBTK_DATA_PATH not set"

```bash
# Configurar manualmente
export GTDBTK_DATA_PATH="/data/database/gtdbtk_251103/gtdbtk_data_release226"
# O agregar al ~/.bashrc
```

---

## 📝 Checklist Final

- [x] Paquete extraído en `/mnt/Part4/Laboratory/pipeline_modular_completo/`
- [x] Directorios creados: `data/`, `output/`, `tmp/`
- [x] Rutas configuradas en `metagenomics_pipeline.sh`
- [x] Temporales configurados: `source setup_tmp_part4.sh`
- [x] Datos copiados a `data/MergedFastq/`
- [x] Espacio verificado: >500 GB disponibles
- [x] Ambientes instalados y verificados
- [x] Variables de entorno configuradas
- [x] Permisos verificados

---

**¡Listo para ejecutar!** 🚀

```bash
cd /mnt/Part4/Laboratory/pipeline_modular_completo
bash metagenomics_pipeline.sh
```

---

**Última actualización:** 3 de diciembre de 2025  
**Ruta configurada:** `/mnt/Part4/Laboratory/tmp`
