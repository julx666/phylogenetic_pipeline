#!/usr/bin/env python3
"""
Script 3: Identifies 1-to-1 orthologous clusters from MMseqs2 clustering results
"""

import pandas as pd
import os
from collections import Counter, defaultdict
import argparse


def extract_genome(gene_name):
    """Extract genome name from gene identifier"""
    parts = gene_name.split('_')
    if len(parts) >= 3:
        return f"{parts[0]}_{parts[1]}"  # e.g. Salmonella_enterica
    elif len(parts) == 2:
        return gene_name  # e.g. Escherichia_coli
    else:
        return gene_name


def parse_clusters(df):
    """Convert representative->member format to list of clusters"""
    
    clusters = defaultdict(set)
    for _, row in df.iterrows():
        clusters[row['representative']].add(row['member'])
    
    # Convert to DataFrame
    cluster_list = []
    for rep, members in clusters.items():
        cluster_list.append({
            'representative': rep,
            'members': list(members),
            'num_genes': len(members)
        })
    
    clusters_df = pd.DataFrame(cluster_list)
    print(f"Created {len(clusters_df)} clusters")
    
    return clusters_df


def identify_orthologs(clusters_df, min_genomes=26):
    """
    Identify 1-to-1 orthologous clusters
    
    Criteria:
    - At least min_genomes different genomes
    - Maximum 1 gene per genome (1-to-1 orthologs)
    """
    
    ortho_clusters = []
    
    for idx, row in clusters_df.iterrows():
        members = row['members']
        
        # Extract genomes
        genomes = [extract_genome(gene) for gene in members]
        genome_counts = Counter(genomes)
        unique_genomes = set(genomes)

        is_one_to_one = all(count == 1 for count in genome_counts.values())
        has_min_genomes = len(unique_genomes) >= min_genomes
        no_duplicates = len(genomes) == len(unique_genomes)
        
        # Check orthology criteria
        if is_one_to_one and has_min_genomes and no_duplicates:
            ortho_clusters.append({
                'representative': row['representative'],
                'members': members,
                'genomes': list(unique_genomes),
                'num_genomes': len(unique_genomes)
            })
  
        
    return ortho_clusters


def save_results(ortho_clusters, output_prefix):
    """Save ortholog results to multiple files"""
    
    # 1. Orthologs 1-1 table (gene level)
    with open(f"{output_prefix}_1_1.tsv", 'w') as f:
        f.write("family_id\tgenome\tgene_id\n")
        for cluster in ortho_clusters:
            for gene in cluster['members']:
                genome = extract_genome(gene)
                f.write(f"{cluster['representative']}\t{genome}\t{gene}\n")
    
    # 2. Genome families (genome composition per family)
    with open(f"{output_prefix}_families.tsv", 'w') as f:
        f.write("family_id\tgenomes\n")
        for cluster in ortho_clusters:
            genomes = ",".join(sorted(cluster['genomes']))
            f.write(f"{cluster['representative']}\t{genomes}\n")
    
    # 3. Representatives list
    with open(f"{output_prefix}_representatives.txt", 'w') as f:
        for cluster in ortho_clusters:
            f.write(f"{cluster['representative']}\n")
    
    print(f"\nSaved results:")
    print(f"  - {output_prefix}_1_1.tsv")
    print(f"  - {output_prefix}_families.tsv")
    print(f"  - {output_prefix}_representatives.txt")


def save_genes_for_alignment(ortho_clusters, output_file):
    """Save all genes from selected ortholog families for alignment"""
    
    with open(output_file, 'w') as f:
        for cluster in ortho_clusters:
            f.write(f"Family: {cluster['representative']} ({cluster['num_genomes']} species)\n")
            for gene in cluster['members']:
                f.write(f"{gene}\n")
            f.write("\n")
    
    total_genes = sum(len(c['members']) for c in ortho_clusters)
    print(f"Saved {total_genes} genes from {len(ortho_clusters)} families to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Identify 1-to-1 orthologous gene families"
    )
    parser.add_argument(
        "-i", "--input",
        default=os.path.join(os.path.dirname(__file__), "..", "results", "all_clusters_mmseqs.tsv"),
        help="Input TSV file from MMseqs2 (default: all_clusters_mmseqs.tsv)"
    )
    parser.add_argument(
        "-o", "--output-prefix",
        default=os.path.join(os.path.dirname(__file__), "..", "results", "orthologs"),
        help="Output prefix for result files (default: orthologs)"
    )
    parser.add_argument(
        "--min-genomes",
        type=int,
        default=26,
        help="Minimum number of genomes for ortholog group (default: 26)"
    )
    parser.add_argument(
        "--genes-output",
        default=os.path.join(os.path.dirname(__file__), "..", "results", "all_genes_for_alignment.txt"),
        help="Output file for genes to align (default: all_genes_for_alignment.txt)"
    )
    
    args = parser.parse_args()
    
    # Read clustering results
    print("Reading clustering results...")
    df = pd.read_csv(args.input, sep='\t', header=None, names=['representative', 'member'])
    print(f"Loaded {len(df):,} cluster relationships")
    print(f"Unique representatives: {df['representative'].nunique():,}")
    print(f"Unique members: {df['member'].nunique():,}")
    
    # Parse clusters
    print("\nParsing clusters...")
    clusters_df = parse_clusters(df)
    
    # Identify orthologs
    print(f"\nIdentifying 1-to-1 orthologs (min {args.min_genomes} genomes)...")
    ortho_clusters = identify_orthologs(clusters_df, args.min_genomes)
    
    print(f"\nRESULTS:")
    print(f"  1-to-1 orthologs: {len(ortho_clusters)}")
    
    # Save all orthologs
    if ortho_clusters:
        save_results(ortho_clusters, args.output_prefix)
        save_genes_for_alignment(ortho_clusters, args.genes_output)

if __name__ == "__main__":
    main()
