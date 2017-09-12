\name{prune_edges}
\alias{prune_edges}
\title{
De-noises ordered differential topological matrix
}
\description{
Second step of DCD: De-noises the ordered differential topological matrix to generate BDs.
}
\usage{
prune_edges(A, min_size_cluster)
}
\arguments{
  \item{A}{
  Ordered noisy differential topological adjacency matrix
}
  \item{min_size_cluster}{
  Minimum size of a community for it to be considered differential, default value is 7.
}
}
\value{
  De-noised ordered topological adjacency matrix
}
\references{
\insertRef{Rpack:bibtex}{Rdpack}
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}
\examples{

#Denoise the noisy ordered adjacency matrix
library(igraph)
library(WGCNA)
g_A <- sample_grg(200, 0.15, torus=FALSE, coords = FALSE)
A <- get.adjacency(g_A)
g_B = permute(g_A,c(sample(20),21:200))
B <- get.adjacency(g_B)

cosine_sim_A <- TOMsimilarity(as.matrix(A),TOMType = "unsigned",TOMDenom = "min");
cosine_sim_B <- TOMsimilarity(as.matrix(B),TOMType = "unsigned",TOMDenom = "min")
edgelist_A <- get.edgelist(g_A);
edgelist_B <- get.edgelist(g_B);
if (is.null(E(g_A)$weight))
{
  edgelist_A <- cbind(edgelist_A,rep(1,nrow(edgelist_A)));
}else
{
  edgelist_A <- cbind(edgelist_A,E(g_A)$weight);
}
if (is.null(E(g_B)$weight))
{
  edgelist_B <- cbind(edgelist_B,rep(1,nrow(edgelist_B)));
}else
{
  edgelist_B <- cbind(edgelist_B,E(g_B)$weight);
}
edgelist_A <- as.data.frame(edgelist_A);
edgelist_B <- as.data.frame(edgelist_B);

mink <- 7
#Noisy Difference in topological matrices
diff_topological_matrix <- abs(cosine_sim_A-cosine_sim_B);

#Order the nodes in topological graph to create block diagonals
ordered_list <- order_topological_matrix(diff_topological_matrix,mink);
temp_output_adjacency <- diff_topological_matrix[ordered_list,ordered_list];

#Perform the greedy deterministic approach to remove spurious edges and keep significant ones
output_adjacency <- prune_edges(temp_output_adjacency,mink);

}
