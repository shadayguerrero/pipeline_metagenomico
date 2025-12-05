#!/usr/bin/env python3
"""
Script para analizar la distribución de Phyla por Reino
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from biom import load_table
import warnings
warnings.filterwarnings('ignore')

# Configurar estilo
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 10

print("="*70)
print("ANÁLISIS DE REINOS Y PHYLA")
print("="*70)

# Cargar archivo BIOM
biom_table = load_table('cuatroc.biom')

# Convertir tabla BIOM a DataFrame
otu_table = biom_table.to_dataframe()

# Obtener taxonomía completa (sin limpiar prefijos primero)
taxonomy_data_raw = []
for obs_id in biom_table.ids(axis='observation'):
    metadata_obs = biom_table.metadata(obs_id, axis='observation')
    if metadata_obs and 'taxonomy' in metadata_obs:
        tax = metadata_obs['taxonomy']
        taxonomy_data_raw.append(tax)
    else:
        taxonomy_data_raw.append(['Unknown'] * 7)

taxonomy_df_raw = pd.DataFrame(taxonomy_data_raw, 
                               index=biom_table.ids(axis='observation'),
                               columns=['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'])

print("\nReinos únicos encontrados (con prefijos):")
print(taxonomy_df_raw['Kingdom'].unique()[:20])

# Limpiar prefijos
taxonomy_df = taxonomy_df_raw.copy()
for col in taxonomy_df.columns:
    taxonomy_df[col] = taxonomy_df[col].apply(lambda x: x[3:] if len(x) > 3 and '__' in x else x)

print("\nReinos únicos encontrados (limpio):")
reinos_unicos = taxonomy_df['Kingdom'].unique()
print(reinos_unicos)
print(f"\nTotal de reinos: {len(reinos_unicos)}")

# Unir con tabla OTU
otu_with_tax = otu_table.copy()
otu_with_tax = otu_with_tax.join(taxonomy_df)

# Filtrar filas sin asignación
otu_with_tax_filtered = otu_with_tax[
    (otu_with_tax['Kingdom'] != '') & 
    (otu_with_tax['Kingdom'] != 'Unknown') &
    (otu_with_tax['Phylum'] != '') & 
    (otu_with_tax['Phylum'] != 'Unknown')
]

print(f"\nOTUs con taxonomía asignada: {len(otu_with_tax_filtered)}")

# Crear tabla de Reino -> Phylum con abundancias
reino_phylum_abundance = otu_with_tax_filtered.groupby(['Kingdom', 'Phylum'])[otu_table.columns].sum()
reino_phylum_total = reino_phylum_abundance.sum(axis=1).reset_index()
reino_phylum_total.columns = ['Kingdom', 'Phylum', 'Total_Abundance']
reino_phylum_total = reino_phylum_total.sort_values('Total_Abundance', ascending=False)

print("\n" + "="*70)
print("DISTRIBUCIÓN DE PHYLA POR REINO")
print("="*70)

for reino in reino_phylum_total['Kingdom'].unique():
    phyla_reino = reino_phylum_total[reino_phylum_total['Kingdom'] == reino]
    total_reino = phyla_reino['Total_Abundance'].sum()
    n_phyla = len(phyla_reino)
    
    print(f"\n{'='*70}")
    print(f"REINO: {reino}")
    print(f"{'='*70}")
    print(f"Total de Phyla: {n_phyla}")
    print(f"Abundancia total: {total_reino:,.0f} reads ({total_reino/otu_table.sum().sum()*100:.2f}%)")
    print(f"\nTop 10 Phyla más abundantes:")
    print("-" * 70)
    
    for idx, row in phyla_reino.head(10).iterrows():
        pct = row['Total_Abundance'] / total_reino * 100
        print(f"  {row['Phylum']:<30} {row['Total_Abundance']:>12,.0f} reads ({pct:>5.1f}%)")

# Guardar tabla completa
reino_phylum_total.to_csv('figuras_python/reino_phylum_distribucion.csv', index=False)
print(f"\n✓ Tabla guardada: figuras_python/reino_phylum_distribucion.csv")

# Crear visualización
print("\n" + "="*70)
print("GENERANDO VISUALIZACIONES")
print("="*70)

# Seleccionar top 30 combinaciones Reino-Phylum
top_combinations = reino_phylum_total.head(30).copy()
top_combinations['Reino_Phylum'] = top_combinations['Kingdom'] + ' - ' + top_combinations['Phylum']

# Gráfico de barras
fig, ax = plt.subplots(figsize=(14, 10))

colors_map = {
    'Bacteria': '#1f77b4',
    'Archaea': '#ff7f0e', 
    'Eukaryota': '#2ca02c',
    'Viruses': '#d62728'
}

# Asignar colores según el reino
colors = [colors_map.get(k, '#7f7f7f') for k in top_combinations['Kingdom']]

bars = ax.barh(range(len(top_combinations)), top_combinations['Total_Abundance'], color=colors)

ax.set_yticks(range(len(top_combinations)))
ax.set_yticklabels(top_combinations['Reino_Phylum'], fontsize=9)
ax.set_xlabel('Abundancia Total (reads)', fontsize=12, fontweight='bold')
ax.set_title('Top 30 Combinaciones Reino-Phylum por Abundancia', fontsize=14, fontweight='bold')
ax.invert_yaxis()

# Añadir leyenda
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=colors_map.get(k, '#7f7f7f'), label=k) 
                   for k in top_combinations['Kingdom'].unique()]
ax.legend(handles=legend_elements, loc='lower right', title='Reino')

plt.tight_layout()
plt.savefig('figuras_python/reino_phylum_top30.png', dpi=300, bbox_inches='tight')
print("✓ Gráfico guardado: figuras_python/reino_phylum_top30.png")
plt.close()

# Crear gráfico de sunburst (treemap simplificado)
fig, ax = plt.subplots(figsize=(14, 10))

# Agrupar por reino
reino_totals = reino_phylum_total.groupby('Kingdom')['Total_Abundance'].sum().sort_values(ascending=False)

# Crear gráfico de barras apiladas mostrando phyla dentro de cada reino
reinos_to_plot = reino_totals.head(5).index  # Top 5 reinos

data_for_plot = []
for reino in reinos_to_plot:
    phyla_data = reino_phylum_total[reino_phylum_total['Kingdom'] == reino].head(10)
    for _, row in phyla_data.iterrows():
        data_for_plot.append({
            'Reino': reino,
            'Phylum': row['Phylum'],
            'Abundancia': row['Total_Abundance']
        })

df_plot = pd.DataFrame(data_for_plot)

# Crear gráfico de barras apiladas horizontal
reino_order = df_plot.groupby('Reino')['Abundancia'].sum().sort_values(ascending=True).index

fig, ax = plt.subplots(figsize=(14, 8))

bottom = np.zeros(len(reino_order))
phyla_unique = df_plot['Phylum'].unique()
colors_phyla = sns.color_palette("tab20", n_colors=len(phyla_unique))
color_map = dict(zip(phyla_unique, colors_phyla))

for phylum in phyla_unique:
    phylum_data = []
    for reino in reino_order:
        value = df_plot[(df_plot['Reino'] == reino) & (df_plot['Phylum'] == phylum)]['Abundancia'].sum()
        phylum_data.append(value)
    
    ax.barh(range(len(reino_order)), phylum_data, left=bottom, label=phylum, color=color_map[phylum])
    bottom += phylum_data

ax.set_yticks(range(len(reino_order)))
ax.set_yticklabels(reino_order, fontsize=12, fontweight='bold')
ax.set_xlabel('Abundancia Total (reads)', fontsize=12, fontweight='bold')
ax.set_title('Composición de Phyla dentro de cada Reino (Top 10 Phyla por Reino)', fontsize=14, fontweight='bold')
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, ncol=2)

plt.tight_layout()
plt.savefig('figuras_python/reino_composicion_phyla.png', dpi=300, bbox_inches='tight')
print("✓ Gráfico guardado: figuras_python/reino_composicion_phyla.png")
plt.close()

# Resumen final
print("\n" + "="*70)
print("RESUMEN FINAL")
print("="*70)

total_reads = otu_table.sum().sum()
print(f"\nTotal de reads analizados: {total_reads:,.0f}")

for reino in reino_totals.index:
    abundancia = reino_totals[reino]
    pct = abundancia / total_reads * 100
    n_phyla = len(reino_phylum_total[reino_phylum_total['Kingdom'] == reino])
    print(f"\n{reino}:")
    print(f"  - Abundancia: {abundancia:,.0f} reads ({pct:.2f}%)")
    print(f"  - Número de Phyla: {n_phyla}")

print("\n" + "="*70)
print("ANÁLISIS COMPLETADO")
print("="*70)

