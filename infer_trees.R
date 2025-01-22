#!/usr/bin/env Rscript

# Install and load necessary packages

install.packages("argparse")
install.packages("phangorn")
install.packages("TreeDist")

library(argparse)
library(phangorn)
library(TreeDist)

# Parse command-line arguments
parser <- ArgumentParser(description="Construct NJ trees from alignment files.")
parser$add_argument("--alignment_dir", required=TRUE, help="Directory containing alignment files.")
parser$add_argument("--output_dir", required=TRUE, help="Directory to write tree files.")
args <- parser$parse_args()

alignment_dir <- args$alignment_dir
output_dir <- args$output_dir

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

run_nj_tree <- function(file_path, output_path) {
    phy_data <- read.phyDat(file_path, format = "fasta", type = "AA")
    dist_matrix <- dist.ml(phy_data, model = "JTT") # distance matrix using JTT protein substitution model
    nj_tree <- NJ(dist_matrix)
    write.tree(nj_tree, file=output_path)
}

# Process all alignment files
alignment_files <- list.files(alignment_dir, full.names=TRUE, pattern="\\.fasta$")
for (alignment_file in alignment_files) {
    base_name <- tools::file_path_sans_ext(basename(alignment_file))
    output_tree <- file.path(output_dir, paste0(base_name, "_NJ_tree.nwk"))
    run_nj_tree(alignment_file, output_tree)
    cat("NJ Tree generated for:", alignment_file, "\n")
}
