\name{calculate_jaccard}
\alias{calculate_jaccard}
\title{
Calculate jaccard coefficient w.r.t. to a hub differential node
}
\description{
Calculate jaccard coefficients of all nodes w.r.t. a hub differential node which help in formation of approximate block diagonals
}
\usage{
calculate_jaccard(nodeid, matrix_A, mink)
}
\arguments{
  \item{nodeid}{
  Id of the hub differential node
}
  \item{matrix_A}{
  Noisy differential topological adjacency matrix
}
  \item{mink}{
  Minimum number of neighbours that each node should have in common with hub node to have a non-zero Jaccard co-efficient.
}
}
\value{
Returns a data-frame with fields: jaccard_coefficient, intersection_length and degree. Here jaccard_coefficient between hub node and every node, intersection_length is number of common nodes and degree represents degree of each node in differential topological matrix.
}
\references{
\insertRef{Rpack:bibtex}{Rdpack}
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}
\examples{
library(igraph)
mink<-7;
g_A <- sample_grg(200, 0.15, torus=FALSE, coords = FALSE)
A <- get.adjacency(g_A)
nodeid <- 1;

#Calculate Jaccard coefficient of all nodes w.r.t. nodeid 1 given the adjacency matrix A
output <- calculate_jaccard(nodeid,A,mink)
}
