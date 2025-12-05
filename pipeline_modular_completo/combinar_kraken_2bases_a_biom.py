#!/usr/bin/env python3
"""
Script para combinar resultados de Kraken2 de dos bases de datos y generar archivo BIOM
VERSIÓN 2: Reconoce correctamente Archaea de GTDB

Este script:
1. Lee los reportes de Kraken2 de dos bases de datos (GTDB y ALT)
2. Combina las abundancias evitando duplicados
3. Clasifica correctamente Archaea usando el linaje completo
4. Genera un archivo BIOM compatible con el script de análisis metagenómico

Uso:
    python3 combinar_kraken_a_biom_v2.py <directorio_kraken> <archivo_salida.biom>

Ejemplo:
    python3 combinar_kraken_a_biom_v2.py kraken/ combined_results.biom
"""

import sys
import os
import json
import numpy as np
from collections import defaultdict
import re

def parse_kraken_report_with_lineage(report_file):
    """
    Parsear un archivo de reporte de Kraken2 manteniendo la jerarquía completa
    
    Formato del reporte:
    % reads | num_reads | num_direct | rank | taxid | name
    """
    taxa_data = {}
    lineages = {}
    stack = []  # Pila para mantener la jerarquía
    
    with open(report_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            
            pct = float(parts[0])
            num_reads = int(parts[1])
            num_direct = int(parts[2])
            rank = parts[3]
            taxid = parts[4]
            name = parts[5]
            
            # Calcular nivel de indentación
            indent_level = (len(name) - len(name.lstrip())) // 2
            name = name.strip()
            
            # Ajustar la pila al nivel actual
            stack = stack[:indent_level]
            
            # Añadir el taxón actual
            current_lineage = stack + [name]
            
            # Guardar lineage
            lineages[taxid] = {
                'name': name,
                'rank': rank,
                'lineage': current_lineage.copy(),
                'reads': num_reads,
                'reads_direct': num_direct,
                'pct': pct
            }
            
            # Guardar solo si tiene reads directos
            if num_direct > 0:
                taxa_data[taxid] = {
                    'name': name,
                    'rank': rank,
                    'reads': num_direct,
                    'pct': pct,
                    'lineage': current_lineage.copy()
                }
            
            stack.append(name)
    
    return taxa_data, lineages

def classify_kingdom_from_lineage(lineage):
    """
    Clasificar el Reino basándose en el linaje completo
    
    Esto es importante para GTDB donde Archaea puede no estar explícito en el nombre del Phylum
    """
    lineage_str = ' -> '.join(lineage)
    
    # Buscar "Archaea" en el linaje
    if 'Archaea' in lineage_str:
        return 'Archaea'
    
    # Buscar "Bacteria" en el linaje
    if 'Bacteria' in lineage_str:
        return 'Bacteria'
    
    # Buscar "Eukaryota" en el linaje
    if 'Eukaryota' in lineage_str or 'Viridiplantae' in lineage_str or 'Metazoa' in lineage_str:
        return 'Eukaryota'
    
    # Buscar indicadores de virus
    if 'Viruses' in lineage_str or 'viricota' in lineage_str.lower():
        return 'Viruses'
    
    # Si no se puede determinar
    return 'Unclassified'

def get_phylum_from_lineage(lineage):
    """
    Extraer el Phylum del linaje
    """
    # Buscar el primer elemento que parezca un Phylum
    # Generalmente es el tercer nivel después de root y dominio
    if len(lineage) >= 3:
        return lineage[2]
    elif len(lineage) >= 2:
        return lineage[1]
    elif len(lineage) >= 1:
        return lineage[0]
    return 'Unknown'

def combine_kraken_reports(gtdb_file, alt_file):
    """
    Combinar dos reportes de Kraken2
    
    Prioridad: GTDB (bacterias/archaea) tiene prioridad, ALT (eucariotas/virus) complementa
    """
    print(f"Leyendo reporte GTDB: {gtdb_file}")
    gtdb_data, gtdb_lineages = parse_kraken_report_with_lineage(gtdb_file)
    
    print(f"Leyendo reporte ALT: {alt_file}")
    alt_data, alt_lineages = parse_kraken_report_with_lineage(alt_file)
    
    # Combinar datos
    combined_data = {}
    
    # Estadísticas
    stats = {
        'Bacteria': 0,
        'Archaea': 0,
        'Eukaryota': 0,
        'Viruses': 0,
        'Unclassified': 0
    }
    
    # Añadir todos los datos de GTDB
    for taxid, data in gtdb_data.items():
        kingdom = classify_kingdom_from_lineage(data['lineage'])
        phylum = get_phylum_from_lineage(data['lineage'])
        
        combined_data[f"GTDB_{taxid}"] = {
            'reads': data['reads'],
            'name': data['name'],
            'rank': data['rank'],
            'lineage': data['lineage'],
            'kingdom': kingdom,
            'phylum': phylum
        }
        
        stats[kingdom] += 1
    
    # Añadir datos de ALT (solo eucariotas, virus, etc.)
    for taxid, data in alt_data.items():
        kingdom = classify_kingdom_from_lineage(data['lineage'])
        phylum = get_phylum_from_lineage(data['lineage'])
        
        # Filtrar solo eucariotas, virus, etc. (no bacterias/archaea)
        if kingdom not in ['Bacteria', 'Archaea']:
            combined_data[f"ALT_{taxid}"] = {
                'reads': data['reads'],
                'name': data['name'],
                'rank': data['rank'],
                'lineage': data['lineage'],
                'kingdom': kingdom,
                'phylum': phylum
            }
            
            stats[kingdom] += 1
    
    print(f"\nTotal de taxones combinados: {len(combined_data)}")
    print(f"  - Bacteria: {stats['Bacteria']}")
    print(f"  - Archaea: {stats['Archaea']}")
    print(f"  - Eukaryota: {stats['Eukaryota']}")
    print(f"  - Viruses: {stats['Viruses']}")
    print(f"  - Unclassified: {stats['Unclassified']}")
    
    return combined_data

def create_biom_from_combined_data(combined_samples, output_file):
    """
    Crear archivo BIOM desde los datos combinados de múltiples muestras
    """
    print("\n" + "="*70)
    print("CREANDO ARCHIVO BIOM")
    print("="*70)
    
    # Recopilar todos los OTUs únicos
    all_otus = set()
    for sample_data in combined_samples.values():
        all_otus.update(sample_data.keys())
    
    all_otus = sorted(list(all_otus))
    sample_ids = sorted(list(combined_samples.keys()))
    
    print(f"\nTotal de OTUs únicos: {len(all_otus)}")
    print(f"Total de muestras: {len(sample_ids)}")
    
    # Estadísticas por reino
    kingdom_counts = {'Bacteria': 0, 'Archaea': 0, 'Eukaryota': 0, 'Viruses': 0, 'Unclassified': 0}
    
    # Crear matriz de abundancias
    matrix = []
    otu_metadata = []
    
    for otu_id in all_otus:
        row = []
        for sample_id in sample_ids:
            if otu_id in combined_samples[sample_id]:
                row.append(combined_samples[sample_id][otu_id]['reads'])
            else:
                row.append(0)
        matrix.append(row)
        
        # Obtener metadata del OTU (de la primera muestra que lo tenga)
        otu_info = None
        for sample_data in combined_samples.values():
            if otu_id in sample_data:
                otu_info = sample_data[otu_id]
                break
        
        if otu_info:
            # Construir taxonomía de 7 niveles
            lineage = otu_info['lineage']
            kingdom = otu_info['kingdom']
            phylum = otu_info['phylum']
            
            # Contar por reino
            kingdom_counts[kingdom] += 1
            
            taxonomy = [''] * 7
            
            # Asignar Reino
            taxonomy[0] = f"k__{kingdom}"
            
            # Asignar Phylum
            if len(lineage) >= 3:
                taxonomy[1] = f"p__{lineage[2]}"
            else:
                taxonomy[1] = f"p__{phylum}"
            
            # Asignar niveles restantes desde el linaje
            for i in range(3, min(len(lineage), 7)):
                prefix = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__'][i]
                taxonomy[i] = f"{prefix}{lineage[i]}"
            
            otu_metadata.append({
                'taxonomy': taxonomy
            })
        else:
            otu_metadata.append({
                'taxonomy': ['k__Unknown'] + [''] * 6
            })
    
    print(f"\nDistribución de OTUs por Reino:")
    for kingdom, count in sorted(kingdom_counts.items(), key=lambda x: x[1], reverse=True):
        pct = (count / len(all_otus)) * 100
        print(f"  - {kingdom}: {count:,} OTUs ({pct:.2f}%)")
    
    # Crear estructura BIOM
    biom_data = {
        "id": "combined_kraken_results",
        "format": "Biological Observation Matrix 1.0.0",
        "format_url": "http://biom-format.org",
        "type": "OTU table",
        "generated_by": "combinar_kraken_a_biom_v2.py",
        "date": "2025-11-03",
        "rows": [
            {
                "id": otu_id,
                "metadata": otu_metadata[i]
            }
            for i, otu_id in enumerate(all_otus)
        ],
        "columns": [
            {
                "id": sample_id,
                "metadata": None
            }
            for sample_id in sample_ids
        ],
        "matrix_type": "sparse",
        "matrix_element_type": "int",
        "shape": [len(all_otus), len(sample_ids)],
        "data": []
    }
    
    # Convertir matriz a formato sparse
    for i, row in enumerate(matrix):
        for j, value in enumerate(row):
            if value > 0:
                biom_data["data"].append([i, j, value])
    
    # Guardar archivo BIOM (formato JSON)
    print(f"\nGuardando archivo BIOM: {output_file}")
    with open(output_file, 'w') as f:
        json.dump(biom_data, f, indent=2)
    
    print(f"✓ Archivo BIOM creado exitosamente")
    print(f"  - {len(all_otus)} OTUs")
    print(f"  - {len(sample_ids)} muestras")
    print(f"  - {len(biom_data['data'])} valores no-cero")
    
    return biom_data

def main(kraken_dir, output_file):
    """
    Función principal
    """
    print("="*70)
    print("COMBINACIÓN DE RESULTADOS DE KRAKEN2 (v2)")
    print("="*70)
    print(f"\nDirectorio de entrada: {kraken_dir}")
    print(f"Archivo de salida: {output_file}")
    
    # Encontrar todos los archivos de reporte
    gtdb_files = {}
    alt_files = {}
    
    for filename in os.listdir(kraken_dir):
        if filename.endswith('_GTDB.report'):
            sample_id = filename.replace('_GTDB.report', '')
            gtdb_files[sample_id] = os.path.join(kraken_dir, filename)
        elif filename.endswith('_ALT.report'):
            sample_id = filename.replace('_ALT.report', '')
            alt_files[sample_id] = os.path.join(kraken_dir, filename)
    
    print(f"\nMuestras encontradas: {len(gtdb_files)}")
    
    # Verificar que cada muestra tenga ambos reportes
    common_samples = set(gtdb_files.keys()) & set(alt_files.keys())
    print(f"Muestras con ambos reportes: {len(common_samples)}")
    
    if len(common_samples) == 0:
        print("\nError: No se encontraron muestras con ambos reportes (GTDB y ALT)")
        return
    
    # Procesar cada muestra
    combined_samples = {}
    
    for sample_id in sorted(common_samples):
        print(f"\n{'='*70}")
        print(f"Procesando muestra: {sample_id}")
        print(f"{'='*70}")
        
        gtdb_file = gtdb_files[sample_id]
        alt_file = alt_files[sample_id]
        
        combined_data = combine_kraken_reports(gtdb_file, alt_file)
        combined_samples[sample_id] = combined_data
    
    # Crear archivo BIOM
    create_biom_from_combined_data(combined_samples, output_file)
    
    print("\n" + "="*70)
    print("PROCESO COMPLETADO")
    print("="*70)
    print(f"\nPuedes usar el archivo '{output_file}' con el script de análisis metagenómico:")
    print(f"  python3 analisis_metagenomico_completo.py {output_file} metadata.csv")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Uso: python3 combinar_kraken_a_biom_v2.py <directorio_kraken> <archivo_salida.biom>")
        print("\nEjemplo:")
        print("  python3 combinar_kraken_a_biom_v2.py kraken/ combined_results.biom")
        sys.exit(1)
    
    kraken_dir = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.exists(kraken_dir):
        print(f"Error: No se encuentra el directorio: {kraken_dir}")
        sys.exit(1)
    
    main(kraken_dir, output_file)

