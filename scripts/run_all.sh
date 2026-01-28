#!/usr/bin/env bash
set -euo pipefail

# Run scripts 01-10 with explicit paths
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(cd "${script_dir}/.." && pwd)"

# Helper function to run and time commands
run_timed() {
  local step_name="$1"
  shift
  echo "========================================="
  echo "Starting: $step_name"
  echo "========================================="
  local start=$(date +%s)
  "$@"
  local end=$(date +%s)
  local duration=$((end - start))
  echo "Completed: $step_name (${duration}s)"
  echo ""
}

# Step functions
step_01() {
  run_timed "Step 01: Download proteomes" \
    python3 "${script_dir}/01_download_proteomes.py" \
    -i "${project_dir}/data/chosen_bacteria.txt" \
    -o "${project_dir}/results/all_bacteria_proteomes.fasta"
}

step_02() {
  run_timed "Step 02: Cluster sequences" \
    python3 "${script_dir}/02_cluster_sequences.py" \
    -i "${project_dir}/results/all_bacteria_proteomes.fasta" \
    -o "${project_dir}/results/all_clusters_mmseqs.tsv"
}

step_03() {
  run_timed "Step 03: Identify orthologs" \
    python3 "${script_dir}/03_identify_orthologs.py" \
    -i "${project_dir}/results/all_clusters_mmseqs.tsv" \
    -o "${project_dir}/results/orthologs" \
    --genes-output "${project_dir}/results/all_genes_for_alignment.txt"
}

step_04() {
  run_timed "Step 04: Extract sequences" \
    python3 "${script_dir}/04_extract_sequences.py" \
    -f "${project_dir}/results/all_bacteria_proteomes.fasta" \
    -g "${project_dir}/results/all_genes_for_alignment.txt" \
    -o "${project_dir}/results/family_fastas"
}

step_05() {
  run_timed "Step 05: Align families" \
    python3 "${script_dir}/05_align_families.py" \
    -i "${project_dir}/results/family_fastas" \
    -o "${project_dir}/results/alignments"
}

step_06() {
  run_timed "Step 06: Build gene trees" \
    python3 "${script_dir}/06_families_trees.py" \
    -i "${project_dir}/results/alignments" \
    -o "${project_dir}/results/all_trees.tree"
}


step_07() {
  run_timed "Step 07: Build greedy consensus tree" \
    python3 "${script_dir}/07_build_greedy_consensus_tree.py" \
    -i "${project_dir}/results/all_trees.tree" \
    -o "${project_dir}/trees_nwk/greedy_consensus.nwk"
}

step_08() {
  run_timed "Step 08: Build consensus tree" \
    Rscript "${script_dir}/08_build_consensus_tree.R"
}

step_09() {
  run_timed "Step 09: Build supertree" \
    python3 "${script_dir}/09_build_supertree.py" \
    -i "${project_dir}/results/all_trees.tree" \
    -o "${project_dir}/trees_nwk/supertree.nwk"
}

step_10() {
  run_timed "Step 10: Validate supertree" \
    Rscript "${script_dir}/10_validate_supertree.R"
}

# Parse optional flags
SKIP_STEP_01=false
START_STEP=1
END_STEP=10

while [[ $# -gt 0 ]]; do
  case "$1" in
    --skip-download)
      SKIP_STEP_01=true
      shift
      ;;
    --steps)
      if [[ $# -lt 2 ]]; then
        echo "ERROR: --steps requires a range argument (e.g., 1-9 or 3-5)" >&2
        exit 1
      fi
      RANGE="$2"
      if [[ "$RANGE" =~ ^([0-9]+)-([0-9]+)$ ]]; then
        START_STEP="${BASH_REMATCH[1]}"
        END_STEP="${BASH_REMATCH[2]}"
        if [[ $START_STEP -lt 1 ]] || [[ $END_STEP -gt 10 ]] || [[ $START_STEP -gt $END_STEP ]]; then
          echo "ERROR: Step range must be between 1-10 and start <= end" >&2
          exit 1
        fi
      else
        echo "ERROR: Invalid step range format. Use 'start-end' (e.g., 3-5)" >&2
        exit 1
      fi
      shift 2
      ;;
    --help|-h)
      echo "Usage: $0 [OPTIONS]"
      echo "Options:"
      echo "  --skip-download     Skip step 01 (download) if output already exists"
      echo "  --steps START-END   Run only steps START to END (e.g., --steps 3-5)"
      echo "  --help              Show this help message"
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      echo "Use --help for usage information" >&2
      exit 1
      ;;
  esac
done

TOTAL_START=$(date +%s)
echo "========================================="
echo "Running pipeline steps $START_STEP to $END_STEP"
echo "========================================="
echo ""

# Execute selected steps
for step_num in $(seq $START_STEP $END_STEP); do
  # Skip step 01 if requested and output already exists
  if [[ $step_num -eq 1 ]] && [[ $SKIP_STEP_01 == true ]]; then
    if [[ -f "${project_dir}/results/all_bacteria_proteomes.fasta" ]]; then
      echo "Output file already exists, skipping step 01 (download)"
      continue
    fi
  fi

  step_func="step_$(printf "%02d" $step_num)"
  if declare -f "$step_func" > /dev/null; then
    $step_func
  else
    echo "No function defined for step $step_num"
    exit 1
  fi
done

TOTAL_END=$(date +%s)
TOTAL_DURATION=$((TOTAL_END - TOTAL_START))
echo "========================================="
echo "Pipeline completed successfully!"
echo "Total time: ${TOTAL_DURATION}s ($((TOTAL_DURATION / 60))m $((TOTAL_DURATION % 60))s)"
echo "========================================="
