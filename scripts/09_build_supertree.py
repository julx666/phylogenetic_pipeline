#!/usr/bin/env python3
"""
Script 7: Build a supertree from gene family trees using ASTRAL
"""

import os
import subprocess
import argparse
import sys
from tempfile import NamedTemporaryFile


def build_supertree(input_file, output_file, astral_jar=None, branch_length_mode=2):
    """Build a supertree from a single file containing all gene family trees"""
    
    if not os.path.exists(input_file):
        print(f"ERROR: Input file not found: {input_file}")
        return False

    with open(input_file, 'r') as infile:
        tree_strings = [line.strip() for line in infile if line.strip()]

    if not tree_strings:
        print(f"ERROR: No trees found in {input_file}")
        return False
    
    print(f"Found {len(tree_strings)} gene family trees")
    

    # Assume input trees already have species names only; just copy to temp file
    print("Copying input trees to temp file for ASTRAL input...")
    temp_files_dir = os.path.join(os.path.dirname(__file__), "..", "temporary_files")
    os.makedirs(temp_files_dir, exist_ok=True)
    with NamedTemporaryFile(mode='w', suffix='.trees', delete=False, dir=temp_files_dir) as f:
        for tree_content in tree_strings:
            f.write(tree_content + '\n')
        trees_list_file = f.name
    
    try:
        # Determine ASTRAL command
        if not astral_jar:
            # Try to find ASTRAL JAR in scripts directory
            script_dir = os.path.dirname(os.path.abspath(__file__))
            default_jar = os.path.join(script_dir, "Astral", "astral.5.7.8.jar")
            if os.path.exists(default_jar):
                astral_jar = default_jar
        
        if astral_jar:
            if not os.path.exists(astral_jar):
                print(f"ERROR: ASTRAL JAR file not found at {astral_jar}")
                return False
            cmd = [
                "java",
                "-Xmx8g",  # Allocate 8GB memory to Java
                "-jar",
                astral_jar,
                "-i",
                trees_list_file,
                "-t",
                str(branch_length_mode),
                "-o",
                output_file
            ]
        else:
            # Try to use ASTRAL from PATH
            cmd = [
                "astral",
                "-i",
                trees_list_file,
                "-t",
                str(branch_length_mode),
                "-o",
                output_file
            ]
        
        print(f"Running ASTRAL with {len(tree_strings)} trees")
        
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
        
        if result.stdout:
            print("ASTRAL output:")
            print(result.stdout)
        
        if not os.path.exists(output_file):
            print("ERROR: Supertree file was not created")
            return False
        
        print(f"\nSupertree successfully created (with ASTRAL support values): {output_file}")
        return True
        
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        print("Make sure ASTRAL is installed and available in PATH, or provide --astral-jar path")
        return False
    except subprocess.CalledProcessError as exc:
        print(f"ERROR: ASTRAL failed with return code {exc.returncode}")
        if exc.stderr:
            print(f"stderr: {exc.stderr}")
        if exc.stdout:
            print(f"stdout: {exc.stdout}")
        return False
    finally:
        # Clean up temporary file
        if os.path.exists(trees_list_file):
            os.remove(trees_list_file)


def main():
    parser = argparse.ArgumentParser(
        description="Build a supertree from gene family trees using ASTRAL"
    )
    parser.add_argument(
        "-i", "--input-file",
        default=os.path.join(os.path.dirname(__file__), "..", "results", "all_trees.tree"),
        help="File containing all gene family trees (one per line)"
    )
    parser.add_argument(
        "-o", "--output-file",
        default=os.path.join(os.path.dirname(__file__), "..", "trees_nwk", "supertree.nwk"),
        help="Output file for the supertree"
    )
    parser.add_argument(
        "--astral-jar",
        default=None,
        help="Path to ASTRAL JAR file (if not in PATH)"
    )
    
    args = parser.parse_args()
    
    success = build_supertree(
        args.input_file,
        args.output_file,
        args.astral_jar,
        2
    )
    
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
