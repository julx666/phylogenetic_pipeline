#!/usr/bin/env Rscript
# Build consensus tree with ape

library(TreeDist)
library(ape)
library(png)
library(phangorn)

# Set output directories for visualizations and trees
output_dir <- "/Users/julaswiatkowska/Desktop/genomika_porownawcza/projekt/trees_visualizations"
nwk_dir <- "/Users/julaswiatkowska/Desktop/genomika_porownawcza/projekt/trees_nwk"

# Create directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!dir.exists(nwk_dir)) {
  dir.create(nwk_dir, recursive = TRUE)
}

# Load trees
all_trees <- read.tree("/Users/julaswiatkowska/Desktop/genomika_porownawcza/projekt/results/all_trees.tree")

# Check what has been loaded
cat("Loaded", length(all_trees), "trees\n")


# Unify and check trees
unify_and_check <- function(trees) {
  # Common taxa
  common <- Reduce(intersect, lapply(trees, function(t) t$tip.label))
  # Break if no common taxa
  if (length(common) == 0) {
    cat("No common taxa!\n")
    return(list())
  }
  # Prune trees
  ref <- sort(common)
  result <- lapply(trees, function(t) {
    if (all(ref %in% t$tip.label)) {
      tr <- keep.tip(t, ref)
      tr$tip.label <- ref
      return(tr)
    }
    return(NULL)
  })
  result <- result[!sapply(result, is.null)]
  return(result)
}

final_trees <- unify_and_check(all_trees)


# Build tree with 50% consensus
consensus_correct <- consensus(final_trees, p = 0.5)

# Support values for trees
supports <- prop.clades(consensus_correct, final_trees)

# Replace NaN → NA
supports[is.nan(supports)] <- NA

# Assign to nodes
consensus_correct$node.label <- round(supports, 3)

# Consistency check
stopifnot(length(consensus_correct$node.label) == consensus_correct$Nnode)

# Plot tree
png(file.path(output_dir, "majority_rule_consensus.png"))
par(mar = c(1, 1, 5, 1))
plot(consensus_correct,
     cex = 0.8,
     main = "Majority-rule consensus (p = 0.5)",
     cex.main = 0.9,
     no.margin = FALSE)
nodelabels(text = consensus_correct$node.label,
           cex = 0.7,
           frame = "none")
invisible(dev.off())

# Save tree
write.tree(consensus_correct,
           file = file.path(nwk_dir, "majority_rule_consensus.nwk"))


# Load trees for comparison
reference <- read.tree("/Users/julaswiatkowska/Desktop/genomika_porownawcza/projekt/trees_nwk/reference_tree.nwk")
timetree <- read.tree("/Users/julaswiatkowska/Desktop/genomika_porownawcza/projekt/trees_nwk/timetree.nwk")
consensus <- read.tree("/Users/julaswiatkowska/Desktop/genomika_porownawcza/projekt/trees_nwk/majority_rule_consensus.nwk")

# Transform both trees to unrooted form
reference_unrooted <- unroot(reference)
consensus_unrooted <- unroot(consensus)
timetree_unrooted <- unroot(timetree)

# Common taxa
common_taxa <- intersect(reference_unrooted$tip.label, consensus_unrooted$tip.label)

# Prune both trees
reference_pruned <- keep.tip(reference_unrooted, common_taxa)
consensus_pruned <- keep.tip(consensus_unrooted, common_taxa)

# RF distance between consensus and reference 
rf_distance <- RF.dist(reference_pruned, consensus_pruned, normalize = FALSE)
rf_normalized <- RF.dist(reference_pruned, consensus_pruned, normalize = TRUE)

cat("Robinson-Foulds distance (Consensus vs Reference):", rf_distance, "\n")
cat("Normalized RF distance (Consensus vs Reference):", rf_normalized, "\n")
cat("Max possible RF distance for", length(common_taxa), "taxons:", 
    2*(length(common_taxa)-3), "\n")

# RF distance between consensus and timetree
common_taxa_timetree <- intersect(consensus_unrooted$tip.label, timetree_unrooted$tip.label)
consensus_pruned_timetree <- keep.tip(consensus_unrooted, common_taxa_timetree)
timetree_pruned <- keep.tip(timetree_unrooted, common_taxa_timetree)

rf_distance_timetree <- RF.dist(consensus_pruned_timetree, timetree_pruned, normalize = FALSE)
rf_normalized_timetree <- RF.dist(consensus_pruned_timetree, timetree_pruned, normalize = TRUE)

cat("\nRobinson-Foulds distance (Consensus vs Timetree):", rf_distance_timetree, "\n")
cat("Normalized RF distance (Consensus vs Timetree):", rf_normalized_timetree, "\n")
cat("Max possible RF distance for", length(common_taxa_timetree), "taxons:", 
    2*(length(common_taxa_timetree)-3), "\n")


# --- Greedy consensus tree analysis ---
# Load greedy consensus tree
greedy_tree_path <- file.path(nwk_dir, "greedy_consensus.nwk")
if (file.exists(greedy_tree_path)) {
  greedy_consensus <- read.tree(greedy_tree_path)
  greedy_consensus_unrooted <- unroot(greedy_consensus)
  
  # Plot greedy consensus tree
  png(file.path(output_dir, "greedy_consensus.png"))
  par(mar = c(1, 1, 5, 1))
  plot(greedy_consensus_unrooted,
       cex = 0.8,
       main = "Greedy consensus",
       cex.main = 0.9,
       no.margin = FALSE)
  nodelabels(text = greedy_consensus_unrooted$node.label,
             cex = 0.7,
             frame = "none")
  invisible(dev.off())

  # RF distance between greedy consensus and reference
  common_taxa_greedy <- intersect(reference_unrooted$tip.label, greedy_consensus_unrooted$tip.label)
  reference_pruned_greedy <- keep.tip(reference_unrooted, common_taxa_greedy)
  greedy_pruned <- keep.tip(greedy_consensus_unrooted, common_taxa_greedy)
  rf_distance_greedy <- RF.dist(reference_pruned_greedy, greedy_pruned, normalize = FALSE)
  rf_normalized_greedy <- RF.dist(reference_pruned_greedy, greedy_pruned, normalize = TRUE)
  cat("\nRobinson-Foulds distance (Greedy Consensus vs Reference):", rf_distance_greedy, "\n")
  cat("Normalized RF distance (Greedy Consensus vs Reference):", rf_normalized_greedy, "\n")
  cat("Max possible RF distance for", length(common_taxa_greedy), "taxons:", 2*(length(common_taxa_greedy)-3), "\n")

  # RF distance between greedy consensus and timetree
  common_taxa_greedy_timetree <- intersect(greedy_consensus_unrooted$tip.label, timetree_unrooted$tip.label)
  greedy_pruned_timetree <- keep.tip(greedy_consensus_unrooted, common_taxa_greedy_timetree)
  timetree_pruned_greedy <- keep.tip(timetree_unrooted, common_taxa_greedy_timetree)
  rf_distance_greedy_timetree <- RF.dist(greedy_pruned_timetree, timetree_pruned_greedy, normalize = FALSE)
  rf_normalized_greedy_timetree <- RF.dist(greedy_pruned_timetree, timetree_pruned_greedy, normalize = TRUE)
  cat("\nRobinson-Foulds distance (Greedy Consensus vs Timetree):", rf_distance_greedy_timetree, "\n")
  cat("Normalized RF distance (Greedy Consensus vs Timetree):", rf_normalized_greedy_timetree, "\n")
  cat("Max possible RF distance for", length(common_taxa_greedy_timetree), "taxons:", 2*(length(common_taxa_greedy_timetree)-3), "\n")
} else {
  cat("\nGreedy consensus tree file not found:", greedy_tree_path, "\n")
}

# Plot timetree
png(file.path(output_dir, "timetree.png"))
par(mar = c(1, 1, 5, 1))
plot(timetree_unrooted,
     cex = 0.8,
     main = "Timetree",
     cex.main = 0.9,
     no.margin = FALSE)
invisible(dev.off())

cat("\nConsensus trees visualizations saved in:", output_dir, "\n")
cat("Files saved:\n")
cat("- majority_rule_consensus.png\n")
cat("- greedy_consensus.png\n")
cat("- timetree.png\n")
cat("\nConsensus trees (Newick format) saved in:", nwk_dir, "\n")
cat("Files saved:\n")
cat("- majority_rule_consensus.nwk\n")
