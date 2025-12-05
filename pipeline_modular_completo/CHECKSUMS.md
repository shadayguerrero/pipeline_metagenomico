# Verificación de Integridad del Paquete

## Archivos Incluidos

Este paquete contiene **32 archivos** esenciales para el funcionamiento del pipeline metagenómico modular.

### Scripts Principales (2)

1. **metagenomics_pipeline.sh** - Script principal con menú interactivo (10 módulos)
2. **setup_micromamba.sh** - Configuración de micromamba

### Scripts por Módulo (7)

3. **run_trimgalore.sh** - Módulo 1: QC & Trimming
4. **run_host_removal.sh** - Módulo 2: Host Removal
5. **run_megahit.sh** - Módulo 3: Assembly
6. **run_binning_fixed.sh** - Módulo 4: Binning (usado por pipeline principal)
7. **run_kraken2.sh** - Módulo 5: Taxonomía de Reads
8. **run_gtdbtk.sh** - Módulo 6: Taxonomía de Bins
9. **run_prokka.sh** - Módulo 7: Anotación

### Scripts de Procesamiento Kraken2 (3)

10. **combinar_kraken_simple_a_biom.py** - Conversión para modo Simple (1 DB)
11. **combinar_kraken_2bases_a_biom.py** - Conversión para modo Dual (2 DBs)
12. **combinar_kraken_3bases_a_biom_v4.py** - Conversión para modo Triple (3 DBs)

### Scripts de Análisis (2)

13. **analisis_sin_metadatos.py** - Análisis sin archivo de metadatos
14. **analisis_metagenomico_completo.py** - Análisis completo con diversidad alpha/beta

### Archivos de Ambientes YAML (13)

15. **01_trimgalore.yaml** - Ambiente para módulo 1
16. **02_host_removal.yaml** - Ambiente para módulo 2
17. **03_megahit.yaml** - Ambiente para módulo 3
18. **04_binning.yaml** - Ambiente para módulo 4
19. **05_kraken2.yaml** - Ambiente para módulo 5
20. **06_gtdbtk.yaml** - Ambiente para módulo 6
21. **07_prokka.yaml** - Ambiente para módulo 7
22. **08_rgi.yaml** - Ambiente para módulo 8
23. **09_antismash.yaml** - Ambiente para módulo 9
24. **10_analysis.yaml** - Ambiente para módulo 10
25. **qc_assembly_env.yaml** - Ambiente combinado para QC y ensamblaje
26. **binning_env.yaml** - Ambiente para binning
27. **tax_annot_env.yaml** - Ambiente para taxonomía y anotación

### Scripts de Instalación (1)

28. **INSTALACION_AMBIENTES.sh** - Instalación automatizada de ambientes micromamba

### Documentación (4)

29. **README.md** - Documentación principal del paquete
30. **GUIA_RAPIDA.md** - Guía de inicio rápido
31. **RESUMEN_IMPLEMENTACION.md** - Documentación técnica detallada
32. **CHECKSUMS.md** - Este archivo

---

## Verificación de Integridad

Para verificar que todos los archivos están presentes:

```bash
# Extraer el paquete
tar -xzf pipeline_modular_completo.tar.gz
cd pipeline_modular_completo

# Verificar número de archivos (debe mostrar 32)
ls -1 | wc -l

# Listar todos los archivos
ls -1
```

### Archivos Esperados (orden alfabético)

```
01_trimgalore.yaml
02_host_removal.yaml
03_megahit.yaml
04_binning.yaml
05_kraken2.yaml
06_gtdbtk.yaml
07_prokka.yaml
08_rgi.yaml
09_antismash.yaml
10_analysis.yaml
CHECKSUMS.md
GUIA_RAPIDA.md
INSTALACION_AMBIENTES.sh
README.md
RESUMEN_IMPLEMENTACION.md
analisis_metagenomico_completo.py
analisis_sin_metadatos.py
binning_env.yaml
combinar_kraken_2bases_a_biom.py
combinar_kraken_3bases_a_biom_v4.py
combinar_kraken_simple_a_biom.py
metagenomics_pipeline.sh
qc_assembly_env.yaml
run_binning_fixed.sh
run_gtdbtk.sh
run_host_removal.sh
run_kraken2.sh
run_megahit.sh
run_prokka.sh
run_trimgalore.sh
setup_micromamba.sh
tax_annot_env.yaml
```

---

## Verificación de Permisos

Los siguientes archivos deben ser ejecutables:

```bash
# Verificar permisos
ls -lh *.sh

# Hacer ejecutables si es necesario
chmod +x *.sh
```

**Archivos que deben ser ejecutables:**
- `metagenomics_pipeline.sh`
- `setup_micromamba.sh`
- `run_trimgalore.sh`
- `run_host_removal.sh`
- `run_megahit.sh`
- `run_binning_fixed.sh`
- `run_kraken2.sh`
- `run_gtdbtk.sh`
- `run_prokka.sh`
- `INSTALACION_AMBIENTES.sh`

---

## Verificación de Sintaxis

### Scripts Bash

```bash
# Verificar sintaxis de todos los scripts bash
for script in *.sh; do
    echo "Verificando $script..."
    bash -n "$script" && echo "✓ OK" || echo "✗ ERROR"
done
```

### Scripts Python

```bash
# Verificar sintaxis de scripts Python
python3 -m py_compile combinar_kraken_simple_a_biom.py
python3 -m py_compile combinar_kraken_2bases_a_biom.py
python3 -m py_compile combinar_kraken_3bases_a_biom_v4.py
python3 -m py_compile analisis_sin_metadatos.py
python3 -m py_compile analisis_metagenomico_completo.py
```

---

## Verificación de Archivos YAML

```bash
# Verificar que todos los YAMLs están presentes
for i in {01..10}; do
    if [ -f "${i}_*.yaml" ]; then
        echo "✓ Ambiente ${i} presente"
    else
        echo "✗ Ambiente ${i} faltante"
    fi
done
```

---

## Información del Paquete

- **Versión:** 1.0
- **Fecha de creación:** 29 de noviembre de 2025
- **Módulos incluidos:** 10 (completo)
- **Archivos totales:** 32
- **Tamaño del paquete:** ~45 KB (comprimido)
- **Tamaño extraído:** ~270 KB

---

## Estructura de Directorios Recomendada

Después de la instalación, tu estructura debería verse así:

```
/data2/shaday/prueba/
├── pipeline_modular_completo/          # Este paquete
│   ├── metagenomics_pipeline.sh        # ← Ejecutar desde aquí
│   ├── *.sh                            # Scripts auxiliares
│   ├── *.py                            # Scripts Python
│   ├── *.yaml                          # Definiciones de ambientes
│   └── *.md                            # Documentación
│
├── MergedFastq/                        # Datos de entrada
│   ├── SRR5936076_1.fastq.gz
│   ├── SRR5936076_2.fastq.gz
│   ├── SRR5936077_1.fastq.gz
│   └── SRR5936077_2.fastq.gz
│
└── output/                             # Resultados (se crea automáticamente)
    ├── trim/
    ├── host_removed/
    ├── megahit_assemblies/
    ├── binning/
    ├── kraken2_*/
    ├── gtdbtk/
    ├── prokka/
    ├── rgi/
    ├── antismash/
    └── analysis/
```

---

## Compatibilidad

- **Sistema Operativo:** Linux (Ubuntu 20.04+, CentOS 7+)
- **Shell:** Bash 4.0+
- **Python:** 3.8+
- **Gestor de ambientes:** Micromamba/Conda

---

## Notas de Versión

### v1.0 (29/11/2025)

**Contenido completo:**
- ✅ 10 scripts de shell (pipeline + módulos individuales)
- ✅ 5 scripts Python (Kraken2 + análisis)
- ✅ 13 archivos YAML (definiciones de ambientes)
- ✅ 4 documentos de ayuda (README, guías, resumen)

**Nuevas características:**
- ✨ Módulo 6: GTDB-Tk (Taxonomía de bins)
- ✨ Módulo 7: Prokka (Anotación funcional)
- ✨ Módulo 8: RGI (Resistencia a antibióticos)
- ✨ Módulo 9: AntiSMASH (Metabolitos secundarios)
- 📚 Documentación completa
- 🛠️ Scripts individuales por módulo
- 🔧 Archivos YAML para instalación de ambientes

---

## Soporte

Para problemas o preguntas:
1. Consultar `README.md` para información general
2. Revisar `GUIA_RAPIDA.md` para casos de uso comunes
3. Leer `RESUMEN_IMPLEMENTACION.md` para detalles técnicos
4. Verificar que todos los 32 archivos estén presentes
5. Comprobar permisos de ejecución en scripts `.sh`

---

**Última actualización:** 29 de noviembre de 2025  
**Archivos verificados:** 32/32 ✅
