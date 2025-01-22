# Install missing packages
list.of.packages <- c('TreeDist', 'ape', 'argparse')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, 'Package'])]
if (length(new.packages)) install.packages(new.packages)

# Load libraries
library(TreeDist)
library(ape)
library(argparse)

# Set up argument parser
parser <- ArgumentParser(description = "Compare trees using shared and unique phylogenetic information.")
parser$add_argument("--input_tree", required = TRUE, help = "File path to tree you want to compare with the reference.")
parser$add_argument("--reference_tree", default = "reference_data/pruned_phylogeny.tre", help = "File path with the reference tree. Default: reference_data/pruned_phylogeny.tre")
args <- parser$parse_args()

# Parse input trees
tree <- read.tree(args$input_tree)
reference_tree <- read.tree(args$reference_tree)

# Change labels to match
# timetree$tip.label <- as.integer(as.factor(timetree$tip.label))
# tree_16s$tip.label <- as.integer(as.factor(tree_16s$tip.label))

# Rotate tree to match the reference
tree <- rotateConstr(tree, reference_tree$tip.label)

# Visualize the matching
VisualizeMatching(
  SharedPhylogeneticInfo, 
  reference_tree, 
  tree, 
  Plot = TreeDistPlot, 
  matchZeros = FALSE
)

# Calculate and print shared and unique information
shared_info <- SharedPhylogeneticInfo(reference_tree, tree)
reference_info <- SplitwiseInfo(reference_tree)
tree_info <- SplitwiseInfo(tree)

cat("Shared information (normalized by reference):", shared_info / reference_info, "\n")
cat("Shared information (normalized by tree):", shared_info / tree_info, "\n")
