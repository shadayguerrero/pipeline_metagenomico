# Guía de Uso del Pipeline en /Part4

## 🎯 Problema

Tu servidor tiene:
- **`/` (root):** 69.9 GB total, **98% usado** (solo 1.1 GB libre) 🔴
- **`/Part4`:** 3.92 TiB disponible ✅

El pipeline necesita **~320 GB** de espacio temporal, especialmente para MEGAHIT y GTDB-Tk.

---

## ✅ Solución Implementada

El pipeline ahora detecta automáticamente `/Part4` y usa ese directorio para archivos temporales en lugar de `/tmp`.

---

## 🚀 Uso Rápido

### Opción 1: Configuración Automática (Recomendado)

El pipeline detecta `/Part4` automáticamente:

```bash
# 1. Ir al directorio del pipeline
cd /Part4/pipeline_modular_completo

# 2. Ejecutar directamente
bash metagenomics_pipeline.sh
```

El pipeline creará automáticamente `/Part4/tmp/general` y lo usará.

---

### Opción 2: Configuración Manual (Más Control)

Si quieres controlar dónde van los temporales:

```bash
# 1. Configurar directorios temporales
source setup_tmp_part4.sh

# 2. Verificar configuración
echo $TMPDIR
# Debe mostrar: /Part4/tmp/general

# 3. Ejecutar pipeline
bash metagenomics_pipeline.sh
```

---

## 📁 Estructura Recomendada en /Part4

```
/Part4/
├── pipeline_modular_completo/     # Pipeline y scripts
│   ├── metagenomics_pipeline.sh
│   ├── run_*.sh
│   └── setup_tmp_part4.sh
│
├── data/                          # Datos de entrada
│   └── MergedFastq/
│       ├── SRR5936076_R1.fastq.gz
│       └── SRR5936076_R2.fastq.gz
│
├── output/                        # Resultados del pipeline
│   ├── trimmed/
│   ├── host_removed/
│   ├── megahit_assemblies/
│   ├── binning/
│   ├── gtdbtk/
│   └── ...
│
└── tmp/                           # Archivos temporales (SE PUEDE BORRAR)
    ├── general/
    ├── megahit/
    ├── gtdbtk/
    ├── concoct/
    └── antismash/
```

---

## 🔧 Configuración Inicial

### Paso 1: Copiar el Pipeline a /Part4

```bash
# Extraer el paquete en /Part4
cd /Part4
tar -xzf pipeline_modular_completo.tar.gz
cd pipeline_modular_completo
```

### Paso 2: Configurar Rutas en el Pipeline

Edita `metagenomics_pipeline.sh` (líneas 27-32):

```bash
# Directorio base del proyecto
PROJECT_DIR="/Part4/data"           # ← CAMBIAR AQUÍ
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Directorios de entrada/salida
INPUT_DIR="${PROJECT_DIR}/MergedFastq"
OUTPUT_DIR="/Part4/output"          # ← CAMBIAR AQUÍ
```

### Paso 3: Crear Directorios

```bash
mkdir -p /Part4/data/MergedFastq
mkdir -p /Part4/output
mkdir -p /Part4/tmp
```

### Paso 4: Copiar Datos

```bash
# Copiar tus reads a /Part4
cp /ruta/original/*.fastq.gz /Part4/data/MergedFastq/
```

---

## 🧪 Verificar Configuración

```bash
# 1. Verificar espacio
df -h /Part4
# Debe mostrar ~3.9 TiB disponible

# 2. Verificar permisos
touch /Part4/tmp/test.txt && rm /Part4/tmp/test.txt
echo "✓ Permisos OK"

# 3. Verificar configuración de temporales
source setup_tmp_part4.sh
# Debe mostrar: ✓ Configuración completada
```

---

## 📊 Uso de Espacio Estimado

Para 2 muestras (SRR5936076, SRR5936077):

| Componente | Ubicación | Espacio |
|------------|-----------|---------|
| **Datos de entrada** | `/Part4/data/` | ~10 GB |
| **Resultados** | `/Part4/output/` | ~50 GB |
| **Temporales MEGAHIT** | `/Part4/tmp/megahit/` | ~100 GB |
| **Temporales GTDB-Tk** | `/Part4/tmp/gtdbtk/` | ~150 GB |
| **Otros temporales** | `/Part4/tmp/` | ~70 GB |
| **Total** | | **~380 GB** |

**Conclusión:** Con 3.92 TiB disponibles, tienes espacio más que suficiente.

---

## 🧹 Limpieza Después de Ejecutar

Los archivos temporales se pueden eliminar después de cada ejecución:

```bash
# Limpiar todos los temporales
rm -rf /Part4/tmp/*

# O limpiar solo temporales antiguos (>7 días)
find /Part4/tmp -type f -mtime +7 -delete
```

**Nota:** NO borres `/Part4/output/` ya que contiene los resultados del análisis.

---

## ⚠️ Monitoreo Durante la Ejecución

### Ver Uso de Espacio en Tiempo Real

```bash
# Terminal 1: Ejecutar pipeline
bash metagenomics_pipeline.sh

# Terminal 2: Monitorear espacio
watch -n 10 'df -h /Part4'
```

### Ver Archivos Temporales

```bash
# Ver tamaño de directorios temporales
du -sh /Part4/tmp/*

# Ver archivos más grandes
find /Part4/tmp -type f -exec du -h {} + | sort -rh | head -20
```

---

## 🚨 Solución de Problemas

### Error: "No space left on device"

**Causa:** Aunque /Part4 tiene espacio, algún proceso está usando `/tmp` (root).

**Solución:**

```bash
# Verificar que TMPDIR está configurado
echo $TMPDIR
# Debe mostrar: /Part4/tmp/general

# Si no está configurado:
source setup_tmp_part4.sh

# Verificar espacio en root
df -h /
# Si está lleno, limpiar:
sudo apt clean
sudo journalctl --vacuum-time=7d
```

### Error: "Permission denied" en /Part4/tmp

**Solución:**

```bash
# Dar permisos
sudo chmod -R 1777 /Part4/tmp

# O cambiar propietario
sudo chown -R $USER:$USER /Part4/tmp
```

### MEGAHIT sigue usando /tmp

**Causa:** Variable `MEGAHIT_TMP` no está configurada.

**Solución:**

```bash
# Configurar manualmente
export MEGAHIT_TMP="/Part4/tmp/megahit"
mkdir -p "${MEGAHIT_TMP}"

# Ejecutar pipeline
bash metagenomics_pipeline.sh
```

---

## 📝 Checklist Antes de Ejecutar

- [ ] Pipeline copiado a `/mnt/Part4/Laboratory/pipeline_modular_completo/`
- [ ] Rutas configuradas en `metagenomics_pipeline.sh`
- [ ] Datos copiados a `/mnt/Part4/Laboratory/data/MergedFastq/`
- [ ] Directorios creados: `/mnt/Part4/Laboratory/output/` y `/mnt/Part4/Laboratory/tmp/`
- [ ] Configuración de temporales: `source setup_tmp_part4.sh`
- [ ] Espacio verificado: `df -h /mnt/Part4` (>500 GB disponibles)
- [ ] Ambientes de micromamba instalados
- [ ] Bases de datos configuradas (GTDB-Tk, Kraken2, etc.)

---

## 🎯 Comando Completo de Ejemplo

```bash
# 1. Ir a /mnt/Part4/Laboratory
cd /mnt/Part4/Laboratory

# 2. Extraer pipeline
tar -xzf pipeline_modular_completo.tar.gz
cd pipeline_modular_completo

# 3. Configurar temporales
source setup_tmp_part4.sh

# 4. Editar rutas (si es necesario)
nano metagenomics_pipeline.sh
# Cambiar PROJECT_DIR y OUTPUT_DIR a /mnt/Part4/Laboratory

# 5. Ejecutar
bash metagenomics_pipeline.sh
```

---

## 📚 Archivos Relacionados

- `setup_tmp_part4.sh` - Script de configuración de temporales
- `ANALISIS_ARCHIVOS_TEMPORALES.md` - Análisis detallado de uso de espacio
- `metagenomics_pipeline.sh` - Pipeline principal (ya configurado)

---

**Última actualización:** 3 de diciembre de 2025  
**Problema:** Root (/) al 98% de uso  
**Solución:** Usar /Part4 para datos, resultados y temporales
