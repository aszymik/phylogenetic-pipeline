library(phangorn)

# Set directories
alignment_dir <- "msa"
output_dir <- "trees"

# Function to read alignment, compute distance matrix, and construct NJ tree
run_nj_tree <- function(file_path, output_path) {
  alignment <- read.dna(file_path, format="fasta")
  dist_matrix <- dist.dna(alignment)
  nj_tree <- NJ(dist_matrix)
  
  # Write the tree to file
  write.tree(nj_tree, file=output_path)
}

# Process all alignment files
alignment_files <- list.files(alignment_dir, full.names = TRUE, pattern = "\\.fasta$")
for (alignment_file in alignment_files) {
  output_tree <- file.path(output_dir, paste0(tools::file_path_sans_ext(basename(alignment_file)), "_NJ_tree.nwk"))
  run_nj_tree(alignment_file, output_tree)
  cat(paste("NJ Tree generated for:", alignment_file, "\n"))
}
