#!/usr/bin/env python3
"""
Script 2: Performs all-vs-all comparison and clustering of protein sequences using MMseqs2
"""

import os
import subprocess
import argparse
import glob
import shutil


def remove_mmseqs_db(db_prefix):
    """Remove all files associated with an MMseqs2 database"""
    for file in glob.glob(f"{db_prefix}*"):
        try:
            os.remove(file)
        except OSError:
            pass


def run_mmseqs2_clustering(input_fasta, output_tsv, temp_dir, min_seq_id=0.8, min_coverage=0.8):
    """
    Run MMseqs2 clustering pipeline
    
    Args:
        input_fasta: Input FASTA file with all proteomes
        output_tsv: Output TSV file with clustering results
        temp_dir: Temporary directory for MMseqs2
        min_seq_id: Minimum sequence identity (default: 0.8 = 80%)
        min_coverage: Minimum alignment coverage (default: 0.8 = 80%)
    """
    
    # Create temp directory
    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(os.path.join(temp_dir, "tmp"), exist_ok=True)
    
    # Define database paths
    proteomy_db = os.path.join(temp_dir, "proteomyDB")
    cluster_db = os.path.join(temp_dir, "clusterDB")
    align_db = os.path.join(temp_dir, "alignDB")
    tmp_subdir = os.path.join(temp_dir, "tmp")
    
    # Remove old databases to avoid conflicts
    remove_mmseqs_db(proteomy_db)
    remove_mmseqs_db(align_db)
    remove_mmseqs_db(cluster_db)
    
    # Create MMseqs2 database from input FASTA
    subprocess.run([
        "mmseqs", "createdb",
        input_fasta,
        proteomy_db
    ], check=True)
    
    # Perform all-vs-all search
    subprocess.run([
        "mmseqs", "search",
        proteomy_db, proteomy_db, # the same base = all-vs-all
        align_db, tmp_subdir,
        "-e", "1e-5",
        "--threads", "4"
    ], check=True)
    
    # Cluster sequences
    subprocess.run([
        "mmseqs", "cluster",
        proteomy_db,
        cluster_db, tmp_subdir,
        "--min-seq-id", str(min_seq_id),
        "-c", str(min_coverage),
        "--cluster-mode", "2",
        "--cov-mode", "3"
    ], check=True)
    
    # Create TSV output
    subprocess.run([
        "mmseqs", "createtsv",
        proteomy_db, proteomy_db,
        cluster_db,
        output_tsv
    ], check=True)
    
    print(f"\nClustering results saved to: {output_tsv}")
    
    # Clean up temporary files
    print("Cleaning up temporary files...")
    shutil.rmtree(temp_dir, ignore_errors=True)
    print("Temporary files removed")


def main():
    parser = argparse.ArgumentParser(
        description="Cluster protein sequences using MMseqs2"
    )
    parser.add_argument(
        "-i", "--input",
        default=os.path.join(os.path.dirname(__file__), "..", "results", "all_bacteria_proteomes.fasta"),
        help="Input FASTA file (default: all_bacteria_proteomes.fasta)"
    )
    parser.add_argument(
        "-o", "--output",
        default=os.path.join(os.path.dirname(__file__), "..", "results", "all_clusters_mmseqs.tsv"),
        help="Output TSV file (default: all_clusters_mmseqs.tsv)"
    )
    parser.add_argument(
        "--min-seq-id",
        type=float,
        default=0.8,
        help="Minimum sequence identity (default: 0.8)"
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.8,
        help="Minimum alignment coverage (default: 0.8)"
    )
    
    args = parser.parse_args()
    
    run_mmseqs2_clustering(
        args.input,
        args.output,
        os.path.join(os.path.dirname(__file__), "..", "temporary_files", "mmseqs_temp"),
        args.min_seq_id,
        args.min_coverage
    )


if __name__ == "__main__":
    main()
