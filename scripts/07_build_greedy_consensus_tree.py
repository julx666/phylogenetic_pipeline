#!/usr/bin/env python3
"""
Script 7b: Build a greedy consensus tree using IQ-TREE
"""

import os
import sys
import argparse
import subprocess


def build_greedy_consensus(input_trees, output_tree):
    """Build a greedy consensus tree using IQ-TREE"""
    if not os.path.exists(input_trees):
        print(f"ERROR: Input trees file not found: {input_trees}")
        return False

    nwk_dir = os.path.dirname(output_tree)
    os.makedirs(nwk_dir, exist_ok=True)
    prefix = os.path.splitext(output_tree)[0]

    cmd = [
        "iqtree",
        "-t", input_trees,
        "-con",
        "-pre", prefix
    ]
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print("ERROR: IQ-TREE failed:")
        print(result.stderr)
        return False

    # IQ-TREE writes consensus tree to a file named 'prefix.contree' in Newick format
    consensus_filename = f"{prefix}.contree"
    if os.path.exists(consensus_filename):
        os.rename(consensus_filename, output_tree)
        print(f"Greedy consensus tree written to {output_tree}")
        return True
    else:
        print("ERROR: Consensus tree not generated. Check IQ-TREE output for errors.")
        return False

def main():
    parser = argparse.ArgumentParser(
        description="Build a greedy consensus tree using IQ-TREE"
    )
    parser.add_argument(
        "-i", "--input-trees",
        default=os.path.join(os.path.dirname(__file__), "..", "results", "all_trees.tree"),
        help="File containing all trees (one per line)"
    )
    parser.add_argument(
        "-o", "--output-tree",
        default=os.path.join(os.path.dirname(__file__), "..", "trees_nwk", "greedy_consensus.nwk"),
        help="Output file for the greedy consensus tree"
    )
    args = parser.parse_args()

    success = build_greedy_consensus(
        args.input_trees,
        args.output_tree
    )
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()

