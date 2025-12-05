#!/usr/bin/env python3
"""
Script para combinar resultados de Kraken2 de TRES bases de datos y generar archivo BIOM
VERSIÓN 4: Maneja correctamente los diferentes sistemas de ranks de GTDB y PlusPFP/EuPathDB

Este script:
1. Lee los reportes de Kraken2 de tres bases de datos (GTDB, PlusPFP, EuPathDB)
2. Maneja correctamente los diferentes sistemas de ranks
3. Combina las abundancias evitando duplicados
4. Clasifica correctamente Archaea y Bacteria
5. Extrae Phyla correctamente
6. Genera un archivo BIOM compatible con el script de análisis metagenómico

Uso:
    python3 combinar_kraken_3bases_a_biom_v4.py <directorio_kraken> <archivo_salida.biom>

Ejemplo:
    python3 combinar_kraken_3bases_a_biom_v4.py report/ combined_results.biom
"""

import sys
import os
import json
import numpy as np
from collections import defaultdict
import re

def parse_kraken_report_flexible(report_file):
    """
    Parsear un archivo de reporte de Kraken2 de manera flexible
    
    Maneja tanto GTDB (R1=Dominio, P=Phylum) como PlusPFP/EuPathDB (R2=Eukaryota, K=Reino, P=Phylum)
    """
    taxa_data = {}
    lineage_stack = []  # Pila para mantener la jerarquía
    
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
            lineage_stack = lineage_stack[:indent_level]
            
            # Añadir el taxón actual a la pila
            lineage_stack.append({
                'name': name,
                'rank': rank,
                'indent': indent_level
            })
            
            # Guardar solo si tiene reads directos
            if num_direct > 0:
                # Construir linaje completo
                full_lineage = [item['name'] for item in lineage_stack]
                
                taxa_data[taxid] = {
                    'name': name,
                    'rank': rank,
                    'reads': num_direct,
                    'pct': pct,
                    'lineage': full_lineage.copy(),
                    'lineage_with_ranks': [(item['name'], item['rank']) for item in lineage_stack]
                }
    
    return taxa_data

def extract_taxonomy_from_lineage(lineage_with_ranks):
    """
    Extraer taxonomía de 7 niveles del linaje con ranks
    
    Maneja tanto GTDB como PlusPFP/EuPathDB
    """
    taxonomy = {
        'Kingdom': '',
        'Phylum': '',
        'Class': '',
        'Order': '',
        'Family': '',
        'Genus': '',
        'Species': ''
    }
    
    for name, rank in lineage_with_ranks:
        # Detectar Reino/Dominio
        if rank in ['R1', 'D', 'K'] or name in ['Bacteria', 'Archaea', 'Eukaryota', 'Viruses']:
            if name in ['Bacteria', 'Archaea', 'Eukaryota', 'Viruses']:
                taxonomy['Kingdom'] = name
            elif rank == 'R2' and name == 'Eukaryota':
                taxonomy['Kingdom'] = 'Eukaryota'
            elif rank == 'K':
                # En PlusPFP, K puede ser un sub-reino como Viridiplantae
                # Pero el Reino principal ya debería estar en R2
                if not taxonomy['Kingdom']:
                    taxonomy['Kingdom'] = name
        
        # Phylum
        if rank == 'P' and not taxonomy['Phylum']:
            taxonomy['Phylum'] = name
        
        # Class
        if rank == 'C' and not taxonomy['Class']:
            taxonomy['Class'] = name
        
        # Order
        if rank == 'O' and not taxonomy['Order']:
            taxonomy['Order'] = name
        
        # Family
        if rank == 'F' and not taxonomy['Family']:
            taxonomy['Family'] = name
        
        # Genus
        if rank == 'G' and not taxonomy['Genus']:
            taxonomy['Genus'] = name
        
        # Species
        if rank == 'S' and not taxonomy['Species']:
            taxonomy['Species'] = name
    
    # Si no se detectó Reino, intentar inferir del linaje
    if not taxonomy['Kingdom']:
        lineage_str = ' -> '.join([name for name, _ in lineage_with_ranks])
        if 'Archaea' in lineage_str:
            taxonomy['Kingdom'] = 'Archaea'
        elif 'Bacteria' in lineage_str:
            taxonomy['Kingdom'] = 'Bacteria'
        elif 'Eukaryota' in lineage_str or 'Viridiplantae' in lineage_str or 'Metazoa' in lineage_str:
            taxonomy['Kingdom'] = 'Eukaryota'
        elif 'Virus' in lineage_str or 'Viruses' in lineage_str:
            taxonomy['Kingdom'] = 'Viruses'
    
    return taxonomy

def classify_kingdom(taxonomy):
    """
    Clasificar el Reino
    """
    kingdom = taxonomy.get('Kingdom', '')
    
    if 'Archaea' in kingdom:
        return 'Archaea'
    elif 'Bacteria' in kingdom:
        return 'Bacteria'
    elif 'Eukaryota' in kingdom or kingdom in ['Viridiplantae', 'Metazoa', 'Fungi']:
        return 'Eukaryota'
    elif 'Virus' in kingdom or 'Viruses' in kingdom:
        return 'Viruses'
    else:
        return 'Unclassified'

def combine_kraken_reports_3bases(gtdb_file, pluspfp_file, eupath_file):
    """
    Combinar tres reportes de Kraken2
    """
    combined_data = {}
    stats = {
        'Bacteria': 0,
        'Archaea': 0,
        'Eukaryota': 0,
        'Viruses': 0,
        'Unclassified': 0
    }
    
    # 1. Añadir todos los datos de GTDB
    print(f"Leyendo reporte GTDB: {os.path.basename(gtdb_file)}")
    gtdb_data = parse_kraken_report_flexible(gtdb_file)
    
    for taxid, data in gtdb_data.items():
        taxonomy = extract_taxonomy_from_lineage(data['lineage_with_ranks'])
        kingdom = classify_kingdom(taxonomy)
        
        combined_data[f"GTDB_{taxid}"] = {
            'reads': data['reads'],
            'name': data['name'],
            'rank': data['rank'],
            'taxonomy': taxonomy,
            'kingdom': kingdom,
            'source': 'GTDB'
        }
        
        stats[kingdom] += 1
    
    # 2. Añadir datos de PlusPFP (solo eucariotas, virus, etc.)
    print(f"Leyendo reporte PlusPFP: {os.path.basename(pluspfp_file)}")
    pluspfp_data = parse_kraken_report_flexible(pluspfp_file)
    
    for taxid, data in pluspfp_data.items():
        taxonomy = extract_taxonomy_from_lineage(data['lineage_with_ranks'])
        kingdom = classify_kingdom(taxonomy)
        
        # Filtrar solo eucariotas, virus, etc. (no bacterias/archaea)
        if kingdom not in ['Bacteria', 'Archaea']:
            combined_data[f"PlusPFP_{taxid}"] = {
                'reads': data['reads'],
                'name': data['name'],
                'rank': data['rank'],
                'taxonomy': taxonomy,
                'kingdom': kingdom,
                'source': 'PlusPFP'
            }
            
            stats[kingdom] += 1
    
    # 3. Añadir datos de EuPathDB (solo eucariotas)
    print(f"Leyendo reporte EuPathDB: {os.path.basename(eupath_file)}")
    eupath_data = parse_kraken_report_flexible(eupath_file)
    
    for taxid, data in eupath_data.items():
        taxonomy = extract_taxonomy_from_lineage(data['lineage_with_ranks'])
        kingdom = classify_kingdom(taxonomy)
        
        # Filtrar solo eucariotas (no bacterias/archaea)
        if kingdom not in ['Bacteria', 'Archaea']:
            combined_data[f"EuPathDB_{taxid}"] = {
                'reads': data['reads'],
                'name': data['name'],
                'rank': data['rank'],
                'taxonomy': taxonomy,
                'kingdom': kingdom,
                'source': 'EuPathDB'
            }
            
            stats[kingdom] += 1
    
    print(f"  Total taxones: {len(combined_data)}")
    print(f"    - Bacteria: {stats['Bacteria']}")
    print(f"    - Archaea: {stats['Archaea']}")
    print(f"    - Eukaryota: {stats['Eukaryota']}")
    print(f"    - Viruses: {stats['Viruses']}")
    print(f"    - Unclassified: {stats['Unclassified']}")
    
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
    source_counts = {'GTDB': 0, 'PlusPFP': 0, 'EuPathDB': 0}
    
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
            tax = otu_info['taxonomy']
            kingdom = otu_info['kingdom']
            source = otu_info['source']
            
            # Contar por reino y fuente
            kingdom_counts[kingdom] += 1
            source_counts[source] += 1
            
            # Construir taxonomía de 7 niveles con prefijos
            taxonomy = [
                f"k__{kingdom}",
                f"p__{tax.get('Phylum', '')}",
                f"c__{tax.get('Class', '')}",
                f"o__{tax.get('Order', '')}",
                f"f__{tax.get('Family', '')}",
                f"g__{tax.get('Genus', '')}",
                f"s__{tax.get('Species', '')}"
            ]
            
            otu_metadata.append({'taxonomy': taxonomy})
        else:
            otu_metadata.append({'taxonomy': ['k__Unknown'] + [''] * 6})
    
    print(f"\nDistribución de OTUs por Reino:")
    for kingdom, count in sorted(kingdom_counts.items(), key=lambda x: x[1], reverse=True):
        pct = (count / len(all_otus)) * 100
        print(f"  - {kingdom}: {count:,} OTUs ({pct:.2f}%)")
    
    print(f"\nDistribución de OTUs por Base de Datos:")
    for source, count in sorted(source_counts.items(), key=lambda x: x[1], reverse=True):
        pct = (count / len(all_otus)) * 100
        print(f"  - {source}: {count:,} OTUs ({pct:.2f}%)")
    
    # Crear estructura BIOM
    biom_data = {
        "id": "combined_kraken_3bases_results",
        "format": "Biological Observation Matrix 1.0.0",
        "format_url": "http://biom-format.org",
        "type": "OTU table",
        "generated_by": "combinar_kraken_3bases_a_biom_v4.py",
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
    print("COMBINACIÓN DE RESULTADOS DE KRAKEN2 (3 BASES - v4)")
    print("="*70)
    print(f"\nDirectorio de entrada: {kraken_dir}")
    print(f"Archivo de salida: {output_file}")
    
    # Encontrar todos los archivos de reporte
    gtdb_files = {}
    pluspfp_files = {}
    eupath_files = {}
    
    for filename in os.listdir(kraken_dir):
        if filename.endswith('_GTDB.report'):
            sample_id = filename.replace('_GTDB.report', '')
            gtdb_files[sample_id] = os.path.join(kraken_dir, filename)
        elif filename.endswith('_PlusPFP.report'):
            sample_id = filename.replace('_PlusPFP.report', '')
            pluspfp_files[sample_id] = os.path.join(kraken_dir, filename)
        elif filename.endswith('_EuPathDB.report'):
            sample_id = filename.replace('_EuPathDB.report', '')
            eupath_files[sample_id] = os.path.join(kraken_dir, filename)
    
    print(f"\nMuestras encontradas:")
    print(f"  - GTDB: {len(gtdb_files)}")
    print(f"  - PlusPFP: {len(pluspfp_files)}")
    print(f"  - EuPathDB: {len(eupath_files)}")
    
    # Verificar que cada muestra tenga los tres reportes
    common_samples = set(gtdb_files.keys()) & set(pluspfp_files.keys()) & set(eupath_files.keys())
    print(f"\nMuestras con los 3 reportes: {len(common_samples)}")
    
    if len(common_samples) == 0:
        print("\nError: No se encontraron muestras con los 3 reportes (GTDB, PlusPFP, EuPathDB)")
        return
    
    # Procesar cada muestra
    combined_samples = {}
    
    for sample_id in sorted(common_samples):
        print(f"\n{'='*70}")
        print(f"Procesando muestra: {sample_id}")
        print(f"{'='*70}")
        
        gtdb_file = gtdb_files[sample_id]
        pluspfp_file = pluspfp_files[sample_id]
        eupath_file = eupath_files[sample_id]
        
        combined_data = combine_kraken_reports_3bases(gtdb_file, pluspfp_file, eupath_file)
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
        print("Uso: python3 combinar_kraken_3bases_a_biom_v4.py <directorio_kraken> <archivo_salida.biom>")
        print("\nEjemplo:")
        print("  python3 combinar_kraken_3bases_a_biom_v4.py report/ combined_results.biom")
        sys.exit(1)
    
    kraken_dir = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.exists(kraken_dir):
        print(f"Error: No se encuentra el directorio: {kraken_dir}")
        sys.exit(1)
    
    main(kraken_dir, output_file)

