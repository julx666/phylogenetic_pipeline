#!/usr/bin/env python3
"""
Script 6: Build phylogenetic trees for gene families
"""

import os
import subprocess
import argparse
import glob
from multiprocessing import Pool


def build_ml_tree(args_tuple):
    """Build a single ML tree with FastTree"""

    alignment_file, output_dir = args_tuple
    family_id = os.path.basename(alignment_file).replace("aligned_", "").replace(".fasta", "")
    tree_file = os.path.join(output_dir, f"tree_{family_id}.tree")

    cmd = [
        "FastTree",
        "-gamma",
        "-lg",
        alignment_file
    ]

    try:
        with open(tree_file, "w") as out_f:
            subprocess.run(
                cmd,
                check=True,
                stdout=out_f,
                stderr=subprocess.PIPE,
                text=True
            )
    except FileNotFoundError:
        return (False, alignment_file, "FastTree not found in PATH")
    except subprocess.CalledProcessError as exc:
        stderr_msg = (exc.stderr or "").strip()
        return (False, alignment_file, stderr_msg if stderr_msg else str(exc))
 
    if not os.path.exists(tree_file) or os.path.getsize(tree_file) == 0:
        return (False, alignment_file, "tree file was not created")

    return (True, alignment_file, None)


def combine_trees(output_dir, combined_filename="all_trees.tree"):
    """Combine all individual tree files into one file and delete individual files"""
    tree_files = sorted(glob.glob(os.path.join(output_dir, "tree_*.tree")))
    
    if not tree_files:
        print("WARNING: No tree files found to combine")
        return
    
    combined_path = os.path.join(output_dir, combined_filename)
    
    with open(combined_path, 'w') as outfile:
        for tree_file in tree_files:
            with open(tree_file, 'r') as infile:
                tree_content = infile.read().strip()
                if tree_content:
                    outfile.write(tree_content + '\n')
    
    print(f"\nCombined {len(tree_files)} trees into {combined_path}")
    
    # Delete individual tree files
    for tree_file in tree_files:
        os.remove(tree_file)
    print(f"Removed {len(tree_files)} individual tree files")
    
    return combined_path


def build_all_trees(input_dir, output_dir, num_processes=4):
    os.makedirs(output_dir, exist_ok=True)

    alignment_files = sorted(glob.glob(os.path.join(input_dir, "aligned_*.fasta")))

    if not alignment_files:
        print(f"ERROR: No aligned FASTA files found in {input_dir}")
        return

    print(f"Found {len(alignment_files)} alignments")
    print(f"Running {num_processes} parallel FastTree processes")

    tree_args = [(f, output_dir) for f in alignment_files]

    success_count = 0
    error_count = 0
    errors = []

    with Pool(processes=num_processes) as pool:
        for i, result in enumerate(pool.imap(build_ml_tree, tree_args), 1):
            success, filename, error = result

            if success:
                success_count += 1
            else:
                error_count += 1
                errors.append((filename, error))
                print(f"ERROR in {filename}: {error}")

            if i % 50 == 0:
                print(f"Progress: {i}/{len(alignment_files)}")

    print("\nTree building complete!")
    print(f"Successful: {success_count}")
    print(f"Failed: {error_count}")

    if errors:
        print("\nFailed families (first 10):")
        for filename, error in errors[:10]:
            print(f"  {filename}: {error}")
    
    # Combine all trees into one file
    return output_dir

def main():
        parser = argparse.ArgumentParser(
            description="Build ML trees for gene families using FastTree and combine into one file"
        )
        parser.add_argument(
            "-i", "--input-dir",
            default=os.path.join(os.path.dirname(__file__), "..", "results", "alignments"),
            help="Directory containing aligned FASTA files"
        )
        parser.add_argument(
            "-o", "--output-file",
            default=os.path.join(os.path.dirname(__file__), "..", "results", "all_trees.tree"),
            help="Output file for combined trees"
        )
        parser.add_argument(
            "-p", "--processes",
            type=int,
            default=4,
            help="Number of parallel FastTree processes"
        )

        args = parser.parse_args()
        
        # Create temporary directory for individual trees
        temp_dir = os.path.join(os.path.dirname(args.output_file), ".temp_trees")

        output_dir = build_all_trees(
            args.input_dir,
            temp_dir,
            args.processes
        )
        
        # Combine trees and move to final location
        if output_dir:
            combined_path = combine_trees(output_dir, "all_trees.tree")
            if combined_path:
                # Move combined file to desired location
                final_path = args.output_file
                os.rename(combined_path, final_path)
                print(f"Final output: {final_path}")
                
                # Remove temporary directory
                try:
                    os.rmdir(temp_dir)
                except OSError:
                    pass

if __name__ == "__main__":
    main()