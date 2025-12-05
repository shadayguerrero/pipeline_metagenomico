# Binning Metagenómico Completo

## Descripción General

El módulo de binning ahora ejecuta **4 herramientas** para generar y refinar bins (genomas metagenómicos ensamblados):

1. **MetaBAT2** - Binning basado en cobertura diferencial y composición tetranucleotídica
2. **MaxBin2** - Binning usando marcadores de genes de copia única
3. **CONCOCT** - Binning con fragmentación de contigs y clustering
4. **DAS Tool** - Refinamiento y selección de los mejores bins de los 3 binnners

---

## Flujo de Trabajo

```
Contigs filtrados (≥500 bp)
    ↓
[1] Indexar con Bowtie2
    ↓
[2] Mapear reads → BAM
    ↓
[3] Calcular profundidad/cobertura
    ↓
    ├─→ [4] MetaBAT2 → bins/
    ├─→ [5] MaxBin2  → bins/
    └─→ [6] CONCOCT  → bins/
         ↓
    [7] DAS Tool (refinamiento)
         ↓
    Bins finales (RECOMENDADO)
```

---

## Estructura de Salida

```
output/binning/${SAMPLE}/
├── mapping/
│   ├── filtered_contigs.fa          # Contigs ≥500 bp
│   ├── contigs_index.*              # Índice Bowtie2
│   ├── mapped_sorted.bam            # Reads mapeados
│   └── *.log                        # Logs de mapeo
│
├── metabat2/
│   ├── bin.1.fa                     # Bins de MetaBAT2
│   ├── bin.2.fa
│   ├── ...
│   ├── depth.txt                    # Profundidad de cobertura
│   └── metabat2.log
│
├── maxbin2/
│   ├── bin.001.fasta                # Bins de MaxBin2
│   ├── bin.002.fasta
│   ├── ...
│   ├── abundance.txt                # Abundancia por contig
│   └── maxbin2.log
│
├── concoct/
│   ├── bins/
│   │   ├── 0.fa                     # Bins de CONCOCT
│   │   ├── 1.fa
│   │   └── ...
│   ├── contigs_10K.fa               # Contigs fragmentados
│   ├── coverage_table.tsv           # Tabla de cobertura
│   └── concoct.log
│
└── dastool/
    ├── DASTool_DASTool_bins/        # ← BINS REFINADOS (USAR ESTOS)
    │   ├── bin_1.fa
    │   ├── bin_2.fa
    │   └── ...
    ├── DASTool_summary.tsv          # Resumen de calidad
    ├── DASTool_scores.pdf           # Gráfico de scores
    └── dastool.log
```

---

## Descripción de Cada Binner

### 1. MetaBAT2

**Método:** Clustering basado en:
- Cobertura diferencial entre muestras
- Composición tetranucleotídica (frecuencia de k-mers de 4 bp)

**Ventajas:**
- Rápido y eficiente
- Funciona bien con múltiples muestras
- Bins de alta calidad

**Parámetros:**
- `--minContig 2500` - Contigs mínimos de 2.5 kb
- `--seed 42` - Semilla para reproducibilidad

**Salida:** `metabat2/bin.*.fa`

---

### 2. MaxBin2

**Método:** Clustering basado en:
- Marcadores de genes de copia única (107 marcadores bacterianos, 40 arqueales)
- Cobertura de contigs
- Composición tetranucleotídica

**Ventajas:**
- Excelente para bacterias y arqueas
- Usa información biológica (marcadores)
- Complementario a MetaBAT2

**Parámetros:**
- `-min_contig_length 2500` - Contigs mínimos de 2.5 kb
- `-thread` - Paralelización

**Salida:** `maxbin2/bin.*.fasta`

---

### 3. CONCOCT

**Método:** 
1. Fragmenta contigs en pedazos de 10 kb
2. Calcula cobertura de cada fragmento
3. Clustering usando composición y cobertura
4. Re-ensambla fragmentos en bins

**Ventajas:**
- Maneja bien contigs largos
- Útil para genomas de baja abundancia
- Diferente estrategia que MetaBAT2/MaxBin2

**Parámetros:**
- `-c 10000` - Fragmentos de 10 kb
- `--merge_last` - Fusiona último fragmento

**Salida:** `concoct/bins/*.fa`

---

### 4. DAS Tool (Dereplicator, Aggregator and Scorer)

**Método:**
1. Recibe bins de los 3 binnners
2. Evalúa calidad de cada bin usando:
   - Genes de copia única (SCG - Single Copy Genes)
   - Completitud
   - Contaminación
3. Selecciona el mejor bin para cada genoma
4. Elimina bins duplicados o de baja calidad

**Ventajas:**
- Combina lo mejor de cada binner
- Reduce redundancia
- Mejora calidad general
- **Produce el conjunto final recomendado**

**Parámetros:**
- `--search_engine diamond` - Usa DIAMOND para búsqueda de genes
- `--write_bins` - Escribe bins refinados

**Salida:** `dastool/DASTool_DASTool_bins/*.fa` ← **USAR ESTOS**

---

## Comparación de Binnners

| Característica | MetaBAT2 | MaxBin2 | CONCOCT | DAS Tool |
|----------------|----------|---------|---------|----------|
| **Velocidad** | Rápido | Medio | Lento | Rápido |
| **Precisión** | Alta | Alta | Media | Muy Alta |
| **Completitud** | Media | Alta | Media | Muy Alta |
| **Contaminación** | Baja | Baja | Media | Muy Baja |
| **Bins típicos** | 10-50 | 10-40 | 20-100 | 5-30 |
| **Mejor para** | Genomas abundantes | Bacterias/Arqueas | Genomas raros | Refinamiento |

---

## Interpretación de Resultados

### Número de Bins Esperado

| Tipo de Muestra | MetaBAT2 | MaxBin2 | CONCOCT | DAS Tool |
|-----------------|----------|---------|---------|----------|
| Baja complejidad | 5-15 | 5-10 | 10-30 | 3-8 |
| Complejidad media | 15-40 | 10-30 | 30-80 | 8-20 |
| Alta complejidad | 40-100+ | 30-60 | 80-200 | 20-50 |

### Calidad de Bins

DAS Tool clasifica bins según:

- **High Quality (HQ)**: Completitud >90%, Contaminación <5%
- **Medium Quality (MQ)**: Completitud >50%, Contaminación <10%
- **Low Quality (LQ)**: Resto

**Recomendación:** Usa solo bins HQ y MQ para análisis posteriores.

---

## Uso

### Ejecutar Binning Completo

```bash
# 1. Activar ambiente
micromamba activate binning

# 2. Ejecutar script
cd /data2/shaday/prueba/pipeline_modular_completo
bash run_binning_fixed.sh
```

### Ejecutar desde el Pipeline Principal

```bash
# Ejecutar pipeline interactivo
bash metagenomics_pipeline.sh

# Seleccionar módulo 4 (Binning)
```

---

## Tiempo de Ejecución Estimado

**Hardware:** 40 cores, 128 GB RAM

| Paso | Tiempo (2 muestras) |
|------|---------------------|
| Filtrado + Mapeo | 30-60 min |
| MetaBAT2 | 10-20 min |
| MaxBin2 | 20-40 min |
| CONCOCT | 40-80 min |
| DAS Tool | 10-20 min |
| **Total** | **2-3.5 horas** |

---

## Siguientes Pasos

Después del binning, usa los bins refinados de DAS Tool para:

1. **Módulo 6 (GTDB-Tk)**: Clasificación taxonómica
2. **Módulo 7 (Prokka)**: Anotación funcional
3. **Módulo 8 (RGI)**: Genes de resistencia
4. **Módulo 9 (AntiSMASH)**: Metabolitos secundarios

**Importante:** Configura los módulos 6-9 para usar:
```bash
BINS_DIR="${OUTPUT_DIR}/binning/${SAMPLE}/dastool/DASTool_DASTool_bins/"
```

---

## Solución de Problemas

### Error: "No bins generados"

**Causa:** Contigs muy cortos o baja cobertura

**Solución:**
```bash
# Reducir tamaño mínimo de bins
MIN_BIN_SIZE=1500 bash run_binning_fixed.sh
```

### Error: "MaxBin2 failed"

**Causa:** Falta `pileup.sh` (de BBMap)

**Solución:**
```bash
micromamba install -n binning bbmap
```

### Error: "DAS Tool no encuentra genes"

**Causa:** Falta DIAMOND

**Solución:**
```bash
micromamba install -n binning diamond
```

### Pocos bins en DAS Tool

**Normal:** DAS Tool es conservador y solo selecciona bins de alta calidad. Si obtienes 5-10 bins refinados de 30-50 bins totales, es esperado.

---

## Referencias

- **MetaBAT2:** Kang et al. (2019) PeerJ 7:e7359
- **MaxBin2:** Wu et al. (2016) Bioinformatics 32(4):605-607
- **CONCOCT:** Alneberg et al. (2014) Nature Methods 11:1144-1146
- **DAS Tool:** Sieber et al. (2018) Nature Microbiology 3:836-843

---

**Última actualización:** 1 de diciembre de 2025  
**Versión:** 2.0 (Binning completo con 4 herramientas)
