library(phangorn)
library(TreeDist)

phy_data <- read.phyDat(..., format = "fasta", type = "AA")
dist_matrix <- dist.ml(phy_data, model = "JTT") # distance matrix using JTT protein substitution model
nj_tree <- nj(as.dist(dist_matrix))
write.tree(nj_tree, file = ...)

