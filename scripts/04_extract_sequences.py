#!/usr/bin/env python3
"""
Script 4: Extract sequences and create family FASTA files
"""

import os
import argparse
import re


def extract_species_name(gene_id):
    """Extract species name from gene ID like 'Dickeya_chrysanthemi_WP_226052722.1' or 'Escherichia_coli_YCN68545.1'"""
    match = re.match(r'([A-Z][a-z]+_[a-z]+)_', gene_id)
    if match:
        return match.group(1)
    return gene_id


def extract_sequences_from_fasta(fasta_file, target_genes):
    """
    Extract sequences for specified genes from FASTA file
    
    Args:
        fasta_file: Input FASTA file
        target_genes: List of gene IDs to extract
    
    Returns:
        Dictionary mapping gene IDs to sequences
    """
    
    sequences = {}
    found = set()
    target_set = set(target_genes)
    
    with open(fasta_file, 'r') as f:
        current_id = None
        current_seq = []
        
        for line in f:
            line = line.rstrip()
            
            if line.startswith('>'):
                # Save previous sequence if it was in target_genes
                if current_id in target_set:
                    sequences[current_id] = ''.join(current_seq)
                    found.add(current_id)
                
                # Extract new ID (first word after '>')
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last sequence
        if current_id in target_set:
            sequences[current_id] = ''.join(current_seq)
            found.add(current_id)
    
    # Report missing genes
    missing = target_set - found
    if missing:
        print(f"WARNING: {len(missing)} genes not found (showing first 5):")
        for gene in list(missing)[:5]:
            print(f"   - {gene}")
    
    print(f"Extracted {len(found)}/{len(target_genes)} genes")
    return sequences


def read_genes_file(genes_file):
    """
    Read genes file and organize by families
    
    Returns:
        List of dictionaries with family info and member genes
    """
    
    families = []
    current_family = None
    
    with open(genes_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            if line.startswith("Family:"):
                # Parse family header
                # Format: "Family: rep_gene (N species)"
                parts = line.split()
                rep = parts[1]
                num_species = int(parts[2].strip('('))
                
                current_family = {
                    'representative': rep,
                    'num_species': num_species,
                    'members': []
                }
                families.append(current_family)
                
            elif line and current_family is not None:
                # Gene ID
                current_family['members'].append(line)
    
    return families


def create_family_fasta_files(families, sequences_dict, output_dir="family_fastas"):
    """
    Create separate FASTA files for each family
    
    Args:
        families: List of family dictionaries
        sequences_dict: Dictionary mapping gene IDs to sequences
        output_dir: Output directory for FASTA files
    
    Returns:
        Tuple of (created_files_list, families_without_genes)
    """
    
    os.makedirs(output_dir, exist_ok=True)
    
    created_files = []
    families_without_genes = 0
    
    for i, family in enumerate(families):
        family_id = family['representative']
        fasta_file = os.path.join(output_dir, f"family_{i+1:04d}.fasta")

        genes_written = 0
        with open(fasta_file, 'w') as f:
            # Write header with only species info
            f.write(f"#Family_{i+1}|Species:{family['num_species']}\n")

            # Write sequences for all genes in family
            for gene in family['members']:
                if gene in sequences_dict:
                    species_name = extract_species_name(gene)
                    f.write(f">{species_name}\n")
                    f.write(f"{sequences_dict[gene]}\n")
                    genes_written += 1
                else:
                    print(f"WARNING: Missing sequence for {gene} in family {i+1}")

        if genes_written > 0:
            created_files.append(fasta_file)
        else:
            families_without_genes += 1
            print(f"WARNING: No sequences written for family {i+1}, removing file")
            os.remove(fasta_file)

    return created_files, families_without_genes


def main():
    parser = argparse.ArgumentParser(
        description="Extract sequences and create family FASTA files"
    )
    parser.add_argument(
        "-f", "--fasta",
        default=os.path.join(os.path.dirname(__file__), "..", "results", "all_bacteria_proteomes.fasta"),
        help="Input FASTA file with all proteomes (default: all_bacteria_proteomes.fasta)"
    )
    parser.add_argument(
        "-g", "--genes",
        default=os.path.join(os.path.dirname(__file__), "..", "results", "all_genes_for_alignment.txt"),
        help="Input file with genes to extract (default: all_genes_for_alignment.txt)"
    )
    parser.add_argument(
        "-o", "--output-dir",
        default=os.path.join(os.path.dirname(__file__), "..", "results", "family_fastas"),
        help="Output directory for family FASTA files (default: family_fastas)"
    )
    
    args = parser.parse_args()
    
    # Read genes file
    print("Reading genes file...")
    families = read_genes_file(args.genes)
    print(f"Found {len(families)} families")
    
    # Collect all genes to extract
    all_genes = []
    for family in families:
        all_genes.extend(family['members'])
    print(f"Total genes to extract: {len(all_genes)}")
    
    # Extract sequences
    sequences_dict = extract_sequences_from_fasta(args.fasta, all_genes)
    
    # Create family FASTA files
    created_files, families_without_genes = create_family_fasta_files(families, sequences_dict, args.output_dir)
    
    # Print summary
    print(f"\nSaved {len(created_files)} family FASTA files to: {args.output_dir}")
    if families_without_genes > 0:
        print(f"({families_without_genes} families skipped - no genes found)")
    

if __name__ == "__main__":
    main()
