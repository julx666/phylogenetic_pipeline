
# Phylogenetic Pipeline

This pipeline downloads bacterial proteomes, clusters genes, extracts 1:1 orthologs, aligns them, builds family trees, and combines them into consensus and supertree with validation. 


## Repository layout
- `data/` — inputs (e.g., `chosen_bacteria.txt`).
- `scripts/` — numbered steps `01`–`10`, helper Bash runner `run_all.sh`.
- `results/` — intermediate outputs (proteomes, clusters, ortholog tables, alignments, combined gene trees).
- `trees_nwk/` — final Newick trees (majority and greedy consensus, supertree, reference, timetree).
- `trees_visualizations/` — PNG plots produced in R steps.
- `requirements.txt` — lists required Python packages for the pipeline.
- `report.pdf` - complex project report
- `phylogenetic_pipeline.pdf` - project presentation
---

## Requirements

### 1. Python
- **Python 3.10+**
- Required Python packages (see `requirements.txt`):
   - `biopython`
   - `pandas`

### 2. R
- **R** (tested with R 4.4.2)
- Required R packages:
   - `TreeDist`
   - `ape`
   - `png`
   - `phangorn`

### 3. External Programs (must be installed and on your PATH):
- **mmseqs2** — for clustering protein sequences
- **mafft** — for multiple sequence alignment
- **FastTree** — for phylogenetic tree inference
- **Java 8+** — required for ASTRAL
- **IQTREE** - required for greedy consensus
- **ASTRAL** — species tree inference (JAR provided at `scripts/Astral/astral.5.7.8.jar` or use your own)

---

## Setup Instructions

### Python environment
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Install external programs (macOS example)
```bash
brew install mmseqs2 mafft fasttree openjdk
# For Linux, use your package manager (apt, yum, etc.) or download binaries from official sites.
```

### R packages
In an R session, run:
```r
install.packages(c("TreeDist", "ape", "png", "phangorn"))
```

---

## Running the pipeline
### Bash runner (steps 1–9)
```bash
cd projekt/scripts
./run_all.sh                        # run everything (defaults to 1-9)
./run_all.sh --steps 3-6            # run a subset
./run_all.sh --steps 3-3            # run just one step
./run_all.sh --skip-download        # reuse existing proteomes
```
The runner uses fixed paths under `projekt/`, produces combined gene trees at `results/all_trees.tree`, consensus trees, the ASTRAL supertree, and validation plots.

## Step-by-step summary
1. **Download proteomes** — `scripts/01_download_proteomes.py` 
   - Input: `data/chosen_bacteria.txt` (one species per line).
   - Output: `results/all_bacteria_proteomes.fasta` (skipped with `--skip-download`).
2. **Cluster sequences (MMseqs2)** — `scripts/02_cluster_sequences.py`
   - Input: proteomes FASTA.
   - Output: `results/all_clusters_mmseqs.tsv` (set `--min-seq-id`, `--min-coverage`).
3. **Identify orthologs** — `scripts/03_identify_orthologs.py`
   - Input: clustering TSV.
   - Outputs: `results/orthologs_1_1.tsv`, `results/orthologs_families.tsv`, `results/orthologs_representatives.txt`, `results/all_genes_for_alignment.txt`.
4. **Extract sequences** — `scripts/04_extract_sequences.py`
   - Inputs: proteomes FASTA, gene list.
   - Output: `results/family_fastas/` with one FASTA per family. 
5. **Align families (MAFFT)** — `scripts/05_align_families.py`
   - Input: `results/family_fastas/`.
   - Output: `results/alignments/` with `aligned_*.fasta`.
6. **Build gene trees (FastTree)** — `scripts/06_families_trees.py`
   - Input: alignments.
   - Output: combined gene trees `results/all_trees.tree` (individual tree files are temporary).
7. **Consensus trees (R)** — `scripts/07_build_consensus_tree.R`
   - Input: `results/all_trees.tree` (species names only).
   - Outputs: `trees_nwk/majority_rule_consensus.nwk`, `trees_nwk/extended_majority_rule_consensus.nwk`, PDFs in `trees_visualizations/`. Uses `trees_nwk/reference_tree.nwk` to report Robinson–Foulds distance.
8. **ASTRAL supertree** — `scripts/08_build_supertree.py`
   - Input: `results/all_trees.tree` (species names only).
   - Output: `trees_nwk/supertree.nwk`
     Use `--astral-jar` to point at a custom JAR.
9. **Validate supertree (R)** — `scripts/09_validate_supertree.R`
   - Inputs: `trees_nwk/reference_tree.nwk`, `trees_nwk/supertree.nwk` 
   - Outputs: RF distance stats and PNG plots in `trees_visualizations/`.

## Tips & troubleshooting
- Make sure `mmseqs`, `mafft`, `FastTree`, and `java` are on `PATH` before running. Test with `which mmseqs` etc.
- ASTRAL needs sufficient heap; adjust `-Xmx` in `scripts/08_build_supertree.py` if runs out of memory.
- R scripts assume `trees_nwk/reference_tree.nwk` exists; replace it with your own reference topology if needed.
- Many scripts accept additional arguments (not just input/output); always check the script's help (`-h`/`--help`) for full usage and options.
