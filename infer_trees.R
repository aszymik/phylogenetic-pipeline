#!/usr/bin/env Rscript

# Install and load necessary packages
install.packages("argparse")
install.packages("phangorn")

library(argparse)
library(phangorn)

# Parse arguments
parser <- ArgumentParser(description="Construct NJ trees with or without bootstrapping from alignment files.")
parser$add_argument("--alignment_dir", required=TRUE, help="Directory containing alignment files.")
parser$add_argument("--output_dir", required=TRUE, help="Directory to write tree files.")
parser$add_argument("--use_bootstrap", action="store_true", help="Include bootstrapping in tree construction.")
parser$add_argument("--bootstrap_replicates", type="integer", default=100, help="Number of bootstrap replicates (only if --use_bootstrap is enabled).")
parser$add_argument("--bootstrap_threshold", type="integer", default=70, help="Minimum bootstrap support to retain branches (only if --use_bootstrap is enabled).")
args <- parser$parse_args()

alignment_dir <- args$alignment_dir
output_dir <- args$output_dir
use_bootstrap <- args$use_bootstrap
bootstrap_replicates <- args$bootstrap_replicates
bootstrap_threshold <- args$bootstrap_threshold

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Generate a basic NJ tree
run_nj_tree <- function(file_path, output_path) {
    # Load alignment
    phy_data <- read.phyDat(file_path, format = "fasta", type = "AA")
    dist_matrix <- dist.ml(phy_data, model = "JTT")  # distance matrix using JTT protein substitution model
    nj_tree <- NJ(dist_matrix)
    write.tree(nj_tree, file = output_path)
}

# Generate NJ tree with bootstrapping
run_nj_with_bootstrap <- function(file_path, output_path, replicates, threshold) {
    phy_data <- read.phyDat(file_path, format = "fasta", type = "AA")
    dist_matrix <- dist.ml(phy_data, model = "JTT")

    # Generate initial NJ tree
    nj_tree <- NJ(dist_matrix)

    # Perform bootstrap analysis
    bootstrap_trees <- bootstrap.phyDat(phy_data, FUN = function(x) NJ(dist.ml(x, model = "JTT")), bs = replicates)
    bs_values <- plotBS(nj_tree, bootstrap_trees, type = "none")

    # Remove branches with low bootstrap support
    for (i in seq_along(bs_values)) {
        if (!is.null(bs_values[i]) && bs_values[i] < threshold) {
            nj_tree$edge.length[i] <- NA
        }
    }
    # Save the tree
    write.tree(nj_tree, file = output_path)
}

# Process all alignment files
alignment_files <- list.files(alignment_dir, full.names=TRUE, pattern="\\.fasta$")
for (alignment_file in alignment_files) {
    base_name <- tools::file_path_sans_ext(basename(alignment_file))
    
    # Determine output tree file name
    if (use_bootstrap) {
        output_tree <- file.path(output_dir, paste0(base_name, "_NJ_tree_with_bootstrap.nwk"))
        run_nj_with_bootstrap(alignment_file, output_tree, bootstrap_replicates, bootstrap_threshold)
        cat("NJ Tree with bootstrapping generated for:", alignment_file, "\n")
    } else {
        output_tree <- file.path(output_dir, paste0(base_name, "_NJ_tree.nwk"))
        run_nj_tree(alignment_file, output_tree)
        cat("NJ Tree generated for:", alignment_file, "\n")
    }
}
