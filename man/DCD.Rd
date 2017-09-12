\name{DCD}
\alias{DCD}
\title{
Differential Community Detection for Paired Biological Networks, for e.g. case vs control
}
\description{
The task of identifying differential sub-networks in paired biological networks (A:control,B:case) can be re-phrased as one of finding dense communities in a single noisy differential topological (DT) graph constructed by taking absolute difference between the topological graphs of A and B. In this software, we propose a fast three-stage approach, namely Differential Community Detection (DCD), to identify differential sub-networks as differential communities in a de-noised version of the DT graph. In the first stage, we iteratively re-order the nodes of the DT graph to determine approximate block diagonals present in the DT adjacency matrix using neighbourhood information of the nodes and Jaccard similarity. In the second stage, the ordered DT adjacency matrix is traversed along the diagonal to remove all the edges associated with a node, if that node has no immediate edges within a window. Finally, we apply community detection methods on this de-noised DT graph to discover differential sub-networks as communities.
}
\usage{
DCD(g_A = sample_grg(200, 0.15, torus = FALSE, coords = FALSE),
    g_B = permute(g_A, c(sample(20), 21:200)), method = "Louvain",
    mink = 7, ground_truth = c(rep(1, 20), rep(0, 180)),
    plot_flag = 1, color = "blue", iter = 1, cores = 1)
}
\arguments{
  \item{g_A}{
  Igraph object representing graph A
}
  \item{g_B}{
  Igraph object representing graph B
}
  \item{method}{
  Community detection method can be either: "Louvain", "Infomap" or "Spectral", default "Louvain".
}
  \item{mink}{
  Minimum number of nodes for a community to be considered differential community, default 7.
}
  \item{ground_truth}{
  In case ground truth community information is available i.e. which nodes are part of differential communities and which are not. For example in 2 networks comprising 1000 nodes, if first 100 nodes are differential then ground_truth = c(rep(1,100),rep(0,900)).
}
  \item{plot_flag}{
  Makes a precision recall plot if plot_flag=1 when ground_truth is available, else set plot_flag to 0.
}
  \item{color}{
  Colour for the area under precision recall curve. Can be "blue", "black", "red" etc.
}
  \item{iter}{
  No of iterations of experiment to run, default iter = 1.
}
  \item{cores}{
  No of cores to use to calculate the Jaccard co-efficient in parallel for the hub node during each iteration of the greedy block diagonal identification approach.
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

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{get_ordered_community_output}}, \code{\link{calculate_jaccard}},
\code{\link{order_topological_matrix}}, \code{\link{prune_edges}}
}
\examples{

%%DCD Example 1 Run%%
DCD_output <- DCD()
%%DCD Example 2 %%
\dontrun{
library(igraph)
g_A <- sample_grg(200, 0.15, torus=FALSE, coords = FALSE)
g_B <- permute(g_A,c(sample(20),21:200))
ground_truth <- c(rep(1,20),rep(0,180))
DCD_output2 <- DCD(g_A = g_A, g_B = g_B, method = "Spectral", mink = 7,
                   ground_truth = ground_truth, plot_flag = 1, color = "red",
                   iter = 1, cores = 1 )
  }
}
