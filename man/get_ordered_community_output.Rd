\name{get_ordered_community_output}
\alias{get_ordered_community_output}
\title{
Get differential communities and nodes belonging to each such community
}
\description{
Third step of DCD: Performs community detection on de-noised ordered differential topological graph.
}
\usage{
get_ordered_community_output(g_output, method, output_adjacency, plot_flags = 0,
                             ground_truth = NULL, color = NULL, iter = 0)
}
\arguments{
  \item{g_output}{
  De-noised ordered differential topological graph.
}
  \item{method}{
  Method can be either "Louvain", "Infomap" or "Spectral", by default "Louvain".
}
  \item{output_adjacency}{
  De-noised ordered differential topological adjacency matrix.
}
  \item{plot_flags}{
 Makes a precision recall plot if plot_flag=1 when ground_truth is available, else set plot_flag to 0.
}
  \item{ground_truth}{
    In case ground truth community information is available i.e. which nodes are part of differential communities and which are not. For example in 2 networks comprising 1000 nodes, if first 100 nodes are differential then ground_truth = c(rep(1,100),rep(0,900)).
}
  \item{color}{
   Colour for the area under precision recall curve. Can be "blue", "black", "red" etc.
}
  \item{iter}{
  No of iterations of experiment to run, default iter = 1.
}
}
\value{
Returns a data frame with NodeIds and Predicted_Label columns. The Predicted_Label corresponds to the Differential ClusterId a node belongs to . All nodes with Predicted_Label = 0 are non-differential nodes or are nodes which are not part of any differential community.
}
\references{
\insertRef{Rpack:bibtex}{Rdpack}
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}
\examples{

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

#Order the nodes in topological graph using greedy approach
ordered_list <- order_topological_matrix(diff_topological_matrix,mink);
temp_output_adjacency <- diff_topological_matrix[ordered_list,ordered_list];

output_adjacency <- prune_edges(temp_output_adjacency,mink);
g_output <- graph_from_adjacency_matrix(output_adjacency,mode=c("undirected"),weighted=TRUE);

#Community detection using Louvain method
louvain_cluster_node_rank <- get_ordered_community_output(g_output,"Louvain",
                          output_adjacency);
output <- cbind(as.numeric(louvain_cluster_node_rank$NodeIds),
                louvain_cluster_node_rank$Predicted_Label)
}
