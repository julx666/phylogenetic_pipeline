#!/usr/bin/env python3
"""
Script 1: Download bacterial proteomes from NCBI
"""

import sys
import time
import gzip
import urllib.request
import os
from Bio import Entrez
import argparse


def download_proteomes(input_file, output_file):
    """
    Download proteomes for organisms listed in input file
    
    Args:
        input_file: File containing organism names (one per line)
        output_file: Output FASTA file for all proteomes
    """
    Entrez.email = "j.swiatkows2@student.uw.edu.pl"
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_file)
    os.makedirs(output_dir, exist_ok=True) 
    
    # Read organism list
    with open(input_file, "r") as f:
        organisms = [line.strip() for line in f if line.strip()]
    
    print(f"Found {len(organisms)} organisms to process")
    
    with open(output_file, "w") as out_file:
        for index, current_species in enumerate(organisms, 1):
            print(f"\n[{index}/{len(organisms)}] Processing: {current_species}")
            
            try:
                # Search for complete genome in assembly database
                search_string = f'"{current_species}"[Organism] AND "complete genome"[Assembly Level] AND "latest RefSeq"[Filter]'
                search_handle = Entrez.esearch(db="assembly", term=search_string, retmax=1)
                search_results = Entrez.read(search_handle)
                search_handle.close()
                
                record_ids = search_results["IdList"]
                
                # Fallback to general search if RefSeq not found
                if not record_ids:
                    print(f"  No RefSeq found, trying general search...")
                    general_search = f'"{current_species}"[Organism] AND "complete genome"[Assembly Level]'
                    search_handle = Entrez.esearch(db="assembly", term=general_search, retmax=1)
                    search_results = Entrez.read(search_handle)
                    search_handle.close()
                    record_ids = search_results["IdList"]
                
                if not record_ids:
                    print(f"  ERROR: No complete genome found")
                    continue
                
                # Get assembly information
                assembly_id = record_ids[0]
                fetch_handle = Entrez.esummary(db="assembly", id=assembly_id)
                assembly_info = Entrez.read(fetch_handle)
                fetch_handle.close()
                
                # Get FTP path for proteome
                doc_summary = assembly_info['DocumentSummarySet']['DocumentSummary'][0]
                ftp_path = doc_summary.get('FtpPath_RefSeq', '') or doc_summary.get('FtpPath_GenBank', '')
                if not ftp_path:
                    print(f"ERROR: No FTP path available")
                    continue
                
                # Build proteome URL
                folder_name = ftp_path.split('/')[-1]
                protein_url = f"{ftp_path}/{folder_name}_protein.faa.gz"
                
                # Download and decompress proteome
                print(f" Downloading proteome...")
                temp_gz = os.path.join(output_dir, "temp_proteom.faa.gz")
                urllib.request.urlretrieve(protein_url, temp_gz)
                
                with gzip.open(temp_gz, 'rt') as gz_file:
                    protein_data = gz_file.read()
                
                # Modify headers to include organism name
                lines = protein_data.strip().split('\n')
                current_header = None
                current_sequence = []
                
                for line in lines:
                    if line.startswith('>'):
                        # Write previous sequence
                        if current_header and current_sequence:
                            out_file.write(current_header + '\n')
                            out_file.write(''.join(current_sequence) + '\n')
                        
                        # Create new header with species name
                        original_id = line[1:].split()[0]
                        new_header = f">{current_species.replace(' ', '_')}_{original_id}"
                        current_header = new_header
                        current_sequence = []
                    elif line.strip():
                        current_sequence.append(line.strip())
                
                # Write last sequence
                if current_header and current_sequence:
                    out_file.write(current_header + '\n')
                    out_file.write(''.join(current_sequence) + '\n')
                
                # Cleanup
                os.remove(temp_gz)
                
                num_proteins = protein_data.count('>')
                print(f"Downloaded: {num_proteins} proteins")
                
                # Rate limiting
                time.sleep(1.0)
                
            except Exception as error:
                print(f"ERROR: {str(error)[:100]}")
    
    print(f"\nProteomes saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Download bacterial proteomes from NCBI"
    )
    parser.add_argument(
        "-i", "--input",
        default=os.path.join(os.path.dirname(__file__), "..","chosen_bacteria.txt"),
        help="Input file with organism names (default: chosen_bacteria.txt)"
    )
    parser.add_argument(
        "-o", "--output",
        default=os.path.join(os.path.dirname(__file__), "..", "results", "all_bacteria_proteomes.fasta"),
        help="Output FASTA file (default: all_bacteria_proteomes.fasta)"
    )
    
    args = parser.parse_args()
    download_proteomes(args.input, args.output)

if __name__ == "__main__":
    main()
