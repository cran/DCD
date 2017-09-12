#    This is an implementation of DCD approach for differential community
#    detection in paired biological networks, in form of an R package.
#    Copyright (C) 2017  Raghvendra Mall

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program, see LICENSE.
DCD <- function(g_A = sample_grg(200, 0.15, torus=FALSE, coords = FALSE), g_B = permute(g_A,c(sample(20),21:200)), method = "Louvain", mink = 7, ground_truth=c(rep(1,20),rep(0,180)), plot_flag=1, color = "blue", iter = 1, cores = 1){

  registerDoParallel(cores)
  #Get adjacency matrix
  adjacencyA <- get.adjacency(g_A)
  if (is.null(colnames(adjacencyA)))
  {
    N <- dim(adjacencyA)[1]
    labels <- paste0("N",c(1:N))
    colnames(adjacencyA) <- labels;
    rownames(adjacencyA) <- labels;
  }
  adjacencyB <- get.adjacency(g_B)
  if (is.null(colnames(adjacencyB)))
  {
    N <- dim(adjacencyB)[1]
    labels <- paste0("N",c(1:N))
    colnames(adjacencyB) <- labels;
    rownames(adjacencyB) <- labels;
  }
  #=================================================================================================

  #Create topological networks
  cosine_sim_A <- TOMsimilarity(as.matrix(adjacencyA),TOMType = "unsigned",TOMDenom = "min");
  cosine_sim_B <- TOMsimilarity(as.matrix(adjacencyB),TOMType = "unsigned",TOMDenom = "min")
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

  #==================================================================================================

  mink=mink;
  #Noisy Difference in topological matrices
  diff_topological_matrix <- abs(cosine_sim_A-cosine_sim_B);

  #Order the nodes in topological graph based on similarity and shortest distance to create block diagonals
  ordered_list <- order_topological_matrix(diff_topological_matrix,mink);
  temp_output_adjacency <- diff_topological_matrix[ordered_list,ordered_list];

  #Perform the greedy deterministic approach to remove spurious edges and keep significant ones
  output_adjacency <- prune_edges(temp_output_adjacency,mink);
  g_output <- graph_from_adjacency_matrix(output_adjacency,mode=c("undirected"),weighted=TRUE);

  #=================================================================================================

  #For a community detection technique
  if (plot_flag>0)
  {
    par(mar=c(5,5,2,2),cex.axis=1.4, cex.lab=1.4, cex.main=1.4, cex.sub=1)
  }
  if (method=="Louvain")
  {
    louvain_cluster_node_rank <- get_ordered_community_output(g_output,method,output_adjacency,plot_flag,ground_truth,color,iter);
    output <- cbind(labels[as.numeric(louvain_cluster_node_rank$NodeIds)],louvain_cluster_node_rank$Predicted_Label)
  }else if (method=="Infomap")
  {
    infomap_cluster_node_rank <- get_ordered_community_output(g_output,method,output_adjacency,plot_flag,ground_truth,color,iter)
    output <- cbind(labels[as.numeric(infomap_cluster_node_rank$NodeIds)],infomap_cluster_node_rank$Predicted_Label)
  }else if (method=="Spectral")
  {
    le_cluster_node_rank <- get_ordered_community_output(g_output,method,output_adjacency,plot_flag,ground_truth,color,iter)
    output <- cbind(labels[as.numeric(le_cluster_node_rank$NodeIds)],le_cluster_node_rank$Predicted_Label)
  }
  output <- as.data.frame(output);
  colnames(output) <- c("NodeIds","Predicted_Label")
  output$NodeIds <- as.character(as.vector(output$NodeIds))
  output$Predicted_Label <- as.numeric(as.vector(output$Predicted_Label))
  return(output)
}
