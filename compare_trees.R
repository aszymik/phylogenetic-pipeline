list.of.packages <- c('TreeDist', 'ape', 'argparse')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]
if(length(new.packages)) install.packages(new.packages)

library('TreeDist')
library('ape')
library('argparse')

parser <- ArgumentParser(description="Compare trees using shared and unique phylogenetic information.")
parser$add_argument("--input_tree", required=TRUE, help="File path to tree you want to compare with the reference.")
parser$add_argument("--reference_tree", default="reference_data/timetree.nwk", help="File path with the reference tree. Default: reference_data/timetree.nwk")
args <- parser$parse_args()

tree <- args$input_tree
reference_tree <- args$reference_tree


# tree <- root(tree, outgroup='Catopuma_badia', resolve.root=TRUE)

# Zmieniamy etykiety liści do wejścia do VisualizeMatching
# timetree$tip.label <- as.integer(as.factor(timetree$tip.label))
# tree_16s$tip.label <- as.integer(as.factor(tree_16s$tip.label))

# Zmiana kolejności liści
tree <- rotateConstr(tree, reference_tree$tip.label)

VisualizeMatching(SharedPhylogeneticInfo, reference_tree, tree, Plot = TreeDistPlot, matchZeros = FALSE)

shared_info <- SharedPhylogeneticInfo(reference_tree, tree)
reference_info <- SplitwiseInfo(reference_tree)
tree_info <- SplitwiseInfo(tree)

shared_info / reference_info
shared_info / tree_info
