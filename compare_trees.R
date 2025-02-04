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
parser$add_argument("--reference_tree", default = "reference_data/timetree.newick", help = "File path with the reference tree. Default: reference_data/pruned_phylogeny.tre")
args <- parser$parse_args()

# Parse input trees
tree <- read.tree(args$input_tree)
if (inherits(tree, "multiPhylo")) tree <- tree[[1]]  # if there is a list of trees, select the first one
print(tree)

reference_tree <- read.tree(args$reference_tree)
if (inherits(reference_tree, "multiPhylo")) reference_tree <- reference_tree[[1]]
print(reference_tree)

# Prune trees to include only common tips
common_tips <- intersect(reference_tree$tip.label, tree$tip.label)
print(length(common_tips))

tree <- drop.tip(tree, setdiff(tree$tip.label, common_tips))
reference_tree <- drop.tip(reference_tree, setdiff(reference_tree$tip.label, common_tips))

# Rotate tree to match the reference
tree <- rotateConstr(tree, reference_tree$tip.label)

# Visualize the matching
VisualizeMatching(
  SharedPhylogeneticInfo, 
  reference_tree, 
  tree, 
  matchZeros = FALSE
)

# Calculate and print shared and unique information
shared_info <- SharedPhylogeneticInfo(reference_tree, tree)
reference_info <- SplitwiseInfo(reference_tree)
tree_info <- SplitwiseInfo(tree)

cat("Reference information:", reference_info, "\n")
cat("Calculated tree information:", tree_info, "\n")
cat("Shared information:", shared_info, "\n")

cat("Shared information (normalized by reference):", shared_info / reference_info, "\n")
cat("Shared information (normalized by tree):", shared_info / tree_info, "\n")
