#!/usr/bin/env python3
"""
Script 5: Align gene families using MAFFT
"""

import os
import subprocess
import argparse
import glob
import tempfile
import shutil
from multiprocessing import Pool


def align_single_family(args_tuple):
    """
    Align a single family using MAFFT
    
    Args:
        args_tuple: (fasta_file, output_dir, mafft_args, tmp_dir)
    
    Returns:
        (success, filename, error_message)
    """
    
    fasta_file, output_dir, mafft_args, tmp_dir = args_tuple
    
    # Create output filename
    base_name = os.path.basename(fasta_file).replace('.fasta', '')
    aligned_file = os.path.join(output_dir, f"aligned_{base_name}.fasta")
    
    try:
        # Build MAFFT command
        cmd = ["mafft"] + mafft_args + [fasta_file]
        
        # Run MAFFT
        with open(aligned_file, 'w') as out_f:
            result = subprocess.run(
                cmd,
                stdout=out_f,
                stderr=subprocess.PIPE,
                text=True,
                env={**os.environ, "MAFFT_TMPDIR": tmp_dir}
            )
        
        # Check if output file was created (MAFFT may write to stderr even on success)
        if os.path.exists(aligned_file) and os.path.getsize(aligned_file) > 0:
            return (True, base_name, None)
        else:
            return (False, base_name, "No output file created")
        
    except Exception as e:
        return (False, base_name, str(e)[:200])


def run_mafft_alignments(input_dir, output_dir):
    """
    Run MAFFT alignments on all family FASTA files

    Args:
        input_dir: Directory with family FASTA files
        output_dir: Directory for aligned FASTA files
    """

    os.makedirs(output_dir, exist_ok=True)

    # Dedicated temp dir for MAFFT scratch files (always cleaned after run)
    tmp_dir = tempfile.mkdtemp(prefix="mafft_tmp_")
    
    # Find all FASTA files
    fasta_files = sorted(glob.glob(os.path.join(input_dir, "*.fasta")))
    
    if not fasta_files:
        print(f"ERROR: No FASTA files found in {input_dir}")
        return
    
    print(f"Found {len(fasta_files)} families to align")
    
    # Use MAFFT auto mode
    mafft_args = ["--auto"]
    threads = 4
    print(f"Using MAFFT mode: auto")
    print(f"Threads: {threads}")
    print(f"MAFFT temp dir: {tmp_dir}")
    
    # Prepare arguments for each alignment
    align_args = [(f, output_dir, mafft_args, tmp_dir) for f in fasta_files]
    
    # Run alignments in parallel
    success_count = 0
    error_count = 0
    
    with Pool(processes=threads) as pool:
        for i, result in enumerate(pool.imap(align_single_family, align_args), 1):
            success, filename, error = result
            
            if success:
                success_count += 1
            else:
                error_count += 1
                print(f"ERROR in {filename}: {error}")
            
            # Progress indicator
            if i % 50 == 0:
                print(f"Progress: {i}/{len(fasta_files)}")
    
    print(f"\nAlignment complete!")
    print(f"Successful: {success_count}")
    print(f"Failed: {error_count}")
    print(f"Alignment saved to: {output_dir}")

    shutil.rmtree(tmp_dir, ignore_errors=True)
    print(f"Removed MAFFT temp dir: {tmp_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Align gene families using MAFFT"
    )
    parser.add_argument(
        "-i", "--input-dir",
        default=os.path.join(os.path.dirname(__file__), "..", "results", "family_fastas"),
        help="Input directory with family FASTA files (default: family_fastas)"
    )
    parser.add_argument(
        "-o", "--output-dir",
        default=os.path.join(os.path.dirname(__file__), "..", "results", "alignments"),
        help="Output directory for aligned files (default: alignments)"
    )
    args = parser.parse_args()
    
    run_mafft_alignments(
        args.input_dir,
        args.output_dir
    )


if __name__ == "__main__":
    main()
