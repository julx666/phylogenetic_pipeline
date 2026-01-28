#!/usr/bin/env Rscript
#Visualize and check the Robinson-Foulds distance between reference tree and supertree


library(TreeDist)
library(ape)
library(png)
library(phangorn)

# Extract numeric support values from node labels
# and scale to branch counts based on number of input trees
extract_and_scale_support <- function(labels, num_trees) {
  if (is.null(labels) || is.null(num_trees)) return(NULL)
  vapply(labels, function(x) {
    m <- regexpr("[0-9]+\\.[0-9]*", x)
    if (m[1] == -1) return(NA_real_)
    prop <- as.numeric(substr(x, m, m + attr(m, "match.length") - 1))
    # Convert proportion to count by multiplying by number of trees
    round(prop * num_trees)
  }, numeric(1))
}


# Set output directory for visualizations
output_dir <- "/Users/julaswiatkowska/Desktop/genomika_porownawcza/projekt/trees_visualizations"

# Create directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load supertree
supertree <- read.tree("/Users/julaswiatkowska/Desktop/genomika_porownawcza/projekt/trees_nwk/supertree.nwk")
timetree <- read.tree("/Users/julaswiatkowska/Desktop/genomika_porownawcza/projekt/trees_nwk/timetree.nwk")
reference <- read.tree("/Users/julaswiatkowska/Desktop/genomika_porownawcza/projekt/trees_nwk/reference_tree.nwk")

# Load input trees to count them for scaling support values
input_trees_file <- "/Users/julaswiatkowska/Desktop/genomika_porownawcza/projekt/results/all_trees.tree"
if (file.exists(input_trees_file)) {
  input_trees_lines <- readLines(input_trees_file)
  num_input_trees <- length(input_trees_lines[input_trees_lines != ""])
  cat("Number of input gene family trees:", num_input_trees, "\n")
} else {
  num_input_trees <- NA
  cat("Warning: Could not find input trees file, support values will not be scaled\n")
}

# Find common taxons
common_taxa <- intersect(reference$tip.label, supertree$tip.label)
cat("Number of common taxons:", length(common_taxa), "\n")

# Find missing taxons
missing_in_supertree <- setdiff(reference$tip.label, supertree$tip.label)
missing_in_reference <- setdiff(supertree$tip.label, reference$tip.label)

if (length(missing_in_supertree) > 0) {
  cat("Missing in supertree:", paste(missing_in_supertree, collapse=", "), "\n")
}
if (length(missing_in_reference) > 0) {
  cat("Missing in reference:", paste(missing_in_reference, collapse=", "), "\n")
}

# Keep common taxons
reference <- keep.tip(reference, common_taxa)
supertree <- keep.tip(supertree, common_taxa)

# Unroot trees for 
reference_unrooted <- unroot(reference)
supertree_unrooted <- unroot(supertree)

# RF distance between supertree and reference
rf_distance <- RF.dist(reference_unrooted, supertree_unrooted, normalize = FALSE, rooted=FALSE)
rf_normalized <- RF.dist(reference_unrooted, supertree_unrooted, normalize = TRUE, rooted=FALSE)

cat("Robinson-Foulds distance (Reference vs Supertree):", rf_distance, "\n")
cat("Normalized RF distance (Reference vs Supertree):", rf_normalized, "\n")
cat("Max possible RF distance for", length(common_taxa), "taxons:", 
    2*(length(common_taxa)-3), "\n")

# RF distance between supertree and timetree
common_taxa_timetree <- intersect(supertree_unrooted$tip.label, timetree$tip.label)
supertree_pruned_timetree <- keep.tip(supertree_unrooted, common_taxa_timetree)
timetree_unrooted <- unroot(timetree)
timetree_pruned <- keep.tip(timetree_unrooted, common_taxa_timetree)

rf_distance_timetree <- RF.dist(supertree_pruned_timetree, timetree_pruned, normalize = FALSE, rooted=FALSE)
rf_normalized_timetree <- RF.dist(supertree_pruned_timetree, timetree_pruned, normalize = TRUE, rooted=FALSE)

cat("\nRobinson-Foulds distance (Supertree vs Timetree):", rf_distance_timetree, "\n")
cat("Normalized RF distance (Supertree vs Timetree):", rf_normalized_timetree, "\n")
cat("Max possible RF distance for", length(common_taxa_timetree), "taxons:", 
    2*(length(common_taxa_timetree)-3), "\n")

# Save visualizations
# 1. Reference tree
png(file = file.path(output_dir, "reference_tree.png"), 
    width = 1200, height = 800, res = 150)
par(mar = c(1, 1, 5, 1))
plot(reference, use.edge.length = FALSE, main = "Reference tree", cex = 0.8)
invisible(dev.off())

# 2. Supertree
png(file = file.path(output_dir, "supertree.png"), 
    width = 1200, height = 800, res = 150)
par(mar = c(1, 1, 5, 1))
plot(supertree, use.edge.length = FALSE, main = "ASTRAL supertree", cex = 0.8)
if (!is.na(num_input_trees)) {
  n_support <- extract_and_scale_support(supertree$node.label, num_input_trees)
  if (!is.null(n_support) && any(!is.na(n_support))) {
    nodelabels(text = as.character(n_support), cex = 0.6, frame = "none")
  }
}
invisible(dev.off())

cat("\nTrees visualizations saved in:", output_dir, "\n")
cat("Files saved:\n")
cat("- reference_tree.png\n")
cat("- supertree.png\n")