#!/usr/bin/env python3
"""
Script de Análisis Metagenómico SIN Metadatos
Genera reportes de abundancias por Reino a partir de archivos BIOM
"""

import sys
import os
from biom import load_table
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from io import BytesIO
import base64

def get_kingdom_from_taxonomy(taxonomy):
    """Extrae el Reino de la taxonomía completa"""
    if not taxonomy or taxonomy == 'Unclassified':
        return 'Unclassified'
    
    # El Reino está en el primer nivel
    levels = [x.strip() for x in taxonomy.split(';')]
    if len(levels) > 0:
        kingdom = levels[0].replace('k__', '').strip()
        if kingdom and kingdom != '':
            return kingdom
    
    return 'Unclassified'

def get_phylum_from_taxonomy(taxonomy):
    """Extrae el Phylum de la taxonomía completa"""
    if not taxonomy or taxonomy == 'Unclassified':
        return 'Unclassified'
    
    # El Phylum está en el segundo nivel
    levels = [x.strip() for x in taxonomy.split(';')]
    if len(levels) > 1:
        phylum = levels[1].replace('p__', '').strip()
        if phylum and phylum != '':
            return phylum
    
    return 'Unclassified'

def create_abundance_plots(kingdom_abs, kingdom_rel, title_prefix=""):
    """Crea gráficos de abundancias absolutas y relativas"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
    
    # Colores por reino
    colors = {
        'Bacteria': '#1f77b4',
        'Archaea': '#ff7f0e',
        'Eukaryota': '#2ca02c',
        'Viruses': '#d62728'
    }
    
    # Gráfico de abundancia absoluta
    kingdom_abs.T.plot(kind='bar', stacked=True, ax=ax1, 
                       color=[colors.get(k, '#7f7f7f') for k in kingdom_abs.index],
                       width=0.8)
    ax1.set_title(f'{title_prefix}Abundancia Absoluta por Reino', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Muestra', fontsize=12)
    ax1.set_ylabel('Abundancia (Conteos)', fontsize=12)
    ax1.legend(title='Reino', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax1.grid(axis='y', alpha=0.3)
    plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    # Gráfico de abundancia relativa
    kingdom_rel.T.plot(kind='bar', stacked=True, ax=ax2,
                       color=[colors.get(k, '#7f7f7f') for k in kingdom_rel.index],
                       width=0.8)
    ax2.set_title(f'{title_prefix}Abundancia Relativa por Reino', fontsize=14, fontweight='bold')
    ax2.set_xlabel('Muestra', fontsize=12)
    ax2.set_ylabel('Abundancia Relativa (%)', fontsize=12)
    ax2.legend(title='Reino', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax2.set_ylim(0, 100)
    ax2.grid(axis='y', alpha=0.3)
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    plt.tight_layout()
    
    # Convertir a base64
    buffer = BytesIO()
    plt.savefig(buffer, format='png', dpi=300, bbox_inches='tight')
    buffer.seek(0)
    img_base64 = base64.b64encode(buffer.read()).decode()
    plt.close()
    
    return img_base64

def create_phylum_plots(phylum_data, kingdom, top_n=15):
    """Crea gráficos de abundancias de Phyla para un Reino específico"""
    if phylum_data.empty:
        return None
    
    # Seleccionar top N phyla
    phylum_sums = phylum_data.sum(axis=1).sort_values(ascending=False)
    top_phyla = phylum_sums.head(top_n).index
    phylum_plot = phylum_data.loc[top_phyla]
    
    # Calcular abundancias relativas
    phylum_rel = phylum_plot.div(phylum_plot.sum(axis=0), axis=1) * 100
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
    
    # Colores
    colors = plt.cm.tab20(np.linspace(0, 1, len(top_phyla)))
    
    # Gráfico de abundancia absoluta
    phylum_plot.T.plot(kind='bar', stacked=True, ax=ax1, color=colors, width=0.8)
    ax1.set_title(f'Abundancia Absoluta - {kingdom} (Top {top_n} Phyla)', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Muestra', fontsize=12)
    ax1.set_ylabel('Abundancia (Conteos)', fontsize=12)
    ax1.legend(title='Phylum', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax1.grid(axis='y', alpha=0.3)
    plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    # Gráfico de abundancia relativa
    phylum_rel.T.plot(kind='bar', stacked=True, ax=ax2, color=colors, width=0.8)
    ax2.set_title(f'Abundancia Relativa - {kingdom} (Top {top_n} Phyla)', fontsize=14, fontweight='bold')
    ax2.set_xlabel('Muestra', fontsize=12)
    ax2.set_ylabel('Abundancia Relativa (% de {})'.format(kingdom), fontsize=12)
    ax2.legend(title='Phylum', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax2.set_ylim(0, 100)
    ax2.grid(axis='y', alpha=0.3)
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    plt.tight_layout()
    
    # Convertir a base64
    buffer = BytesIO()
    plt.savefig(buffer, format='png', dpi=300, bbox_inches='tight')
    buffer.seek(0)
    img_base64 = base64.b64encode(buffer.read()).decode()
    plt.close()
    
    return img_base64

def generate_html_report(biom_file, output_dir):
    """Genera reporte HTML sin metadatos"""
    print("Cargando archivo BIOM...")
    table = load_table(biom_file)
    
    # Extraer datos
    sample_ids = table.ids(axis='sample')
    obs_ids = table.ids(axis='observation')
    
    print(f"Muestras: {len(sample_ids)}")
    print(f"OTUs: {len(obs_ids)}")
    
    # Crear DataFrame de taxonomía
    taxonomy_data = []
    for obs_id in obs_ids:
        metadata = table.metadata(obs_id, axis='observation')
        taxonomy = metadata.get('taxonomy', 'Unclassified') if metadata else 'Unclassified'
        taxonomy_data.append({
            'OTU_ID': obs_id,
            'Taxonomy': taxonomy,
            'Kingdom': get_kingdom_from_taxonomy(taxonomy),
            'Phylum': get_phylum_from_taxonomy(taxonomy)
        })
    
    taxonomy_df = pd.DataFrame(taxonomy_data)
    
    # Crear matriz de abundancias
    abundance_matrix = table.matrix_data.toarray()
    abundance_df = pd.DataFrame(abundance_matrix, index=obs_ids, columns=sample_ids)
    
    # Combinar con taxonomía
    abundance_df = abundance_df.join(taxonomy_df.set_index('OTU_ID'))
    
    print("\n1. CALCULANDO ABUNDANCIAS POR REINO")
    
    # Agrupar por Reino
    kingdom_abs = abundance_df.groupby('Kingdom')[sample_ids].sum()
    
    # Filtrar solo reinos clasificados para abundancia relativa
    classified_kingdoms = ['Bacteria', 'Archaea', 'Eukaryota', 'Viruses']
    kingdom_classified = kingdom_abs.loc[kingdom_abs.index.isin(classified_kingdoms)]
    kingdom_rel = kingdom_classified.div(kingdom_classified.sum(axis=0), axis=1) * 100
    
    # Crear gráfico de reinos
    print("   Generando gráficos de Reino...")
    kingdom_plot = create_abundance_plots(kingdom_classified, kingdom_rel)
    
    # Guardar tablas
    kingdom_abs.to_csv(f'{output_dir}/abundancias_absolutas_reino.csv')
    kingdom_rel.to_csv(f'{output_dir}/abundancias_relativas_reino.csv')
    
    print("\n2. CALCULANDO ABUNDANCIAS DE PHYLA POR REINO")
    
    phylum_plots = {}
    for kingdom in classified_kingdoms:
        if kingdom in abundance_df['Kingdom'].values:
            print(f"   Procesando: {kingdom}")
            kingdom_data = abundance_df[abundance_df['Kingdom'] == kingdom]
            phylum_abs = kingdom_data.groupby('Phylum')[sample_ids].sum()
            
            if not phylum_abs.empty:
                plot = create_phylum_plots(phylum_abs, kingdom)
                if plot:
                    phylum_plots[kingdom] = plot
                    
                    # Guardar tablas
                    phylum_abs.to_csv(f'{output_dir}/abundancias_absolutas_{kingdom.lower()}.csv')
                    phylum_rel = phylum_abs.div(phylum_abs.sum(axis=0), axis=1) * 100
                    phylum_rel.to_csv(f'{output_dir}/abundancias_relativas_{kingdom.lower()}.csv')
    
    print("\n3. GENERANDO REPORTE HTML")
    
    # Crear HTML
    html_content = f"""
<!DOCTYPE html>
<html lang="es">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Reporte Metagenómico - Abundancias por Reino</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}
        .header h1 {{
            margin: 0;
            font-size: 2.5em;
        }}
        .header p {{
            margin: 10px 0 0 0;
            opacity: 0.9;
        }}
        .section {{
            background: white;
            padding: 30px;
            margin-bottom: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .section h2 {{
            color: #667eea;
            border-bottom: 3px solid #667eea;
            padding-bottom: 10px;
            margin-top: 0;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
        }}
        .stat-value {{
            font-size: 2.5em;
            font-weight: bold;
            margin: 10px 0;
        }}
        .stat-label {{
            font-size: 0.9em;
            opacity: 0.9;
        }}
        img {{
            max-width: 100%;
            height: auto;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background-color: #667eea;
            color: white;
        }}
        tr:hover {{
            background-color: #f5f5f5;
        }}
        .footer {{
            text-align: center;
            padding: 20px;
            color: #666;
            font-size: 0.9em;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>📊 Reporte Metagenómico</h1>
        <p>Análisis de Abundancias por Reino y Phylum</p>
    </div>
    
    <div class="section">
        <h2>📋 Resumen General</h2>
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-label">Total de Muestras</div>
                <div class="stat-value">{len(sample_ids)}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Total de OTUs</div>
                <div class="stat-value">{len(obs_ids):,}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Reinos Detectados</div>
                <div class="stat-value">{len(kingdom_classified)}</div>
            </div>
        </div>
        
        <h3>Muestras Analizadas</h3>
        <p><strong>{', '.join(sample_ids)}</strong></p>
    </div>
    
    <div class="section">
        <h2>👑 Abundancias por Reino</h2>
        <p>Distribución de organismos en los cuatro reinos principales: Bacteria, Archaea, Eukaryota y Viruses.</p>
        <img src="data:image/png;base64,{kingdom_plot}" alt="Abundancias por Reino">
        
        <h3>Tabla de Abundancias por Reino</h3>
        {kingdom_abs.T.to_html(classes='table', border=0)}
    </div>
"""
    
    # Añadir secciones de Phyla
    for kingdom, plot in phylum_plots.items():
        html_content += f"""
    <div class="section">
        <h2>🧬 Phyla de {kingdom}</h2>
        <p>Composición de Phyla dentro del reino {kingdom}.</p>
        <img src="data:image/png;base64,{plot}" alt="Phyla de {kingdom}">
    </div>
"""
    
    html_content += """
    <div class="footer">
        <p>Reporte generado automáticamente por el Pipeline Metagenómico</p>
        <p>Los archivos CSV con los datos completos están disponibles en el directorio de salida</p>
    </div>
</body>
</html>
"""
    
    # Guardar HTML
    output_file = f'{output_dir}/reporte_abundancias.html'
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"\n✓ Reporte HTML generado: {output_file}")
    print(f"✓ Tablas CSV guardadas en: {output_dir}/")
    
    return output_file

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Uso: python analisis_sin_metadatos.py <archivo.biom> <directorio_salida>")
        sys.exit(1)
    
    biom_file = sys.argv[1]
    output_dir = sys.argv[2]
    
    if not os.path.exists(biom_file):
        print(f"Error: No se encontró el archivo {biom_file}")
        sys.exit(1)
    
    os.makedirs(output_dir, exist_ok=True)
    
    generate_html_report(biom_file, output_dir)

