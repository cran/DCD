get_ordered_community_output <- function(g_output,method,output_adjacency,plot_flags=0,ground_truth=NULL,color=NULL,iter=0)
{
  if (method=="Louvain")
  {
    community_output <- cluster_louvain(g_output)
  }
  else if (method=="Infomap")
  {
    community_output <- cluster_infomap(g_output,E(g_output)$weights,nb.trials=100)
  }
  else if (method=="Spectral")
  {
    community_output <- cluster_leading_eigen(g_output,weights = E(g_output)$weights)
  }
  comm_membership <- membership(community_output)
  cids <- unique(comm_membership)

  #Get volume of all clusters and sorted contribution of each node in a cluster
  cluster_volume <- NULL
  node_ranking <- NULL
  for (i in 1:length(cids))
  {
    cid <- cids[i]
    nodeids <- which(comm_membership==cid)

    #Get contribution of each node
    if (length(nodeids)>1)
    {
      node_contributions <- colSums(output_adjacency[nodeids,nodeids]);

      #Get volume of each cluster
      cluster_contributions <- sum(node_contributions);
    }else{
      node_contributions <- sum(output_adjacency[nodeids,nodeids]);
      cluster_contributions <- node_contributions;
    }
    temp1 <- cbind(cid,cluster_contributions);
    cluster_volume <- rbind(cluster_volume,temp1);

    cluster_contributions <- rep(cluster_contributions,length(node_contributions));
    temp2 <- cbind(nodeids,rep(cid,length(cluster_contributions)));
    node_ranking <- rbind(node_ranking,temp2);
  }

  #Get the node contributions
  node_ranking <- as.data.frame(node_ranking);

  #Get cluster volume
  cluster_volume <- as.data.frame(cluster_volume);

  #Rank clusters based on volume
  final_ranking <- NULL

  #Sort clusters by volume
  cluster_volume <- cluster_volume[order(cluster_volume[,2],decreasing = TRUE),];
  for (i in 1:nrow(cluster_volume))
  {
    cid <- cluster_volume[i,1];
    temp <- node_ranking[node_ranking[,2]==cid,2];
    if (length(temp)>1)
    {
      temp <- rep(i,length(temp))
    }
    else
    {
      temp <- rep(0,length(temp))
    }

    final_ranking <- rbind(final_ranking,cbind(node_ranking[node_ranking[,2]==cid,],temp));
  }

  rownames(final_ranking) <- NULL
  colnames(final_ranking) <- c("NodeIds","ClusterIds","Predicted_Label");
  final_ranking <- as.data.frame(final_ranking);

  ##Sort based on volumne of each cluster
  output_sorted_per_cluster_rank <- final_ranking$Predicted_Label;

  if (plot_flags>0)
  {
    #In order to plot
    pred_cluster <- prediction(output_sorted_per_cluster_rank,ground_truth[final_ranking[,1]]);
    perf_cluster <- performance(pred_cluster,"prec","rec")

    if (iter==1)
    {
      plot(perf_cluster,col=color,lwd=2,lty=3,xlim=c(0,1),ylim=c(0,1.01));
    }
    else if (iter>1)
    {
      plot(perf_cluster,col=color,add=TRUE,lwd=2,lty=3,xlim=c(0,1),ylim=c(0,1.01));
    }
  }
  return(final_ranking)
}

calculate_jaccard <- function(nodeid,matrix_A,mink)
{
  jacc_df <- NULL;
  vec1 <- matrix_A[nodeid,];
  temp_indices <- which(vec1>0)
  nrows <- nrow(matrix_A)
  i <- 0;
  jacc_df <- foreach(i=1:nrows, .combine = rbind, .inorder=TRUE) %dopar%
  {
    temp2_indices <- which(matrix_A[i,]>0);
    intersect_length <- length(intersect(temp_indices,temp2_indices));
    if (intersect_length>=mink)
    {
      union_length <- length(union(temp_indices,temp2_indices));
      value <- intersect_length/union_length
    }else{
      value <- 0;
    }
    output <- cbind(value, intersect_length, length(temp2_indices))
    output
  }
  jacc_df <- as.data.frame(jacc_df);
  colnames(jacc_df) <- c("jaccard_coefficient","intersection_length","degree");
  return(jacc_df);
}

#Function to order the nodes in topological matrix
order_topological_matrix <- function(A,mink)
{
  main_degree_A <- rowSums(A);
  N <- nrow(A);
  ordered_list <- NULL;
  rem_indices <- c(1:N);
  #Convert matrix A to binary matrix
  B = A;
  B[B>0] <- 1;
  main_degree_B <- rowSums(B);
  main_degree_B <- main_degree_B+1;
  rm(B)
  gc();

  #Perform the ordering
  prev_jacc_df <- NULL;
  degree_A <- main_degree_A*log(main_degree_B);
  while(length(ordered_list) < N )
  {
    #Get node with highest degree
    temp_nodeid <- which.max(degree_A);

    #If all the nodes are zero degree nodes then stop
    if (degree_A[temp_nodeid]<=mink)
    {
      ordered_list <- c(ordered_list,rem_indices);
      break;
    }

    nodeid <- rem_indices[temp_nodeid];

    #Calculate jaccard co-efficient w.r.t. the hub
    jacc_df <- calculate_jaccard(nodeid,A,mink);
    indices <- order(jacc_df$jaccard_coefficient,decreasing=T)
    new_jacc_df <- jacc_df[indices,];

    #Get the non-zero jaccard ids
    non_zero_ordered_indices <- indices[new_jacc_df$jaccard_coefficient>0];
    non_zero_ordered_indices <- non_zero_ordered_indices[order(jacc_df[non_zero_ordered_indices,]$degree,decreasing=TRUE)]
    rm(new_jacc_df)
    gc()

    #Jaccard index of hub with itself should be one so >0 indices should be > 1
    if (length(non_zero_ordered_indices)>1)
    {
      if (is.null(prev_jacc_df))
      {
        ordered_list <- c(ordered_list,non_zero_ordered_indices);
      }
      else{
        ids_not_in_common <- setdiff(non_zero_ordered_indices,ordered_list);
        ids_in_common <- intersect(ordered_list,non_zero_ordered_indices);
        ids_in_common_to_consider <- ids_in_common[jacc_df$jaccard_coefficient[ids_in_common]>prev_jacc_df$jaccard_coefficient[ids_in_common]]
        revised_ordered_list <- setdiff(ordered_list,ids_in_common_to_consider);

        if (!is.null(ids_in_common_to_consider))
        {
          ids_in_common_to_consider <- rev(ids_in_common_to_consider);
          revised_ordered_list <- c(revised_ordered_list,ids_in_common_to_consider);
        }
        revised_ordered_list <- c(revised_ordered_list,ids_not_in_common);
        ordered_list <- revised_ordered_list;
      }
    }
    else
    {
      #If no more non-zero index left then break out of loop
      ordered_list <- c(ordered_list,rem_indices);
      break;
    }

    #Left over indices
    rem_indices <- setdiff(c(1:N),ordered_list);

    #If length of remaining indices is 0 then break out of the while loop
    if (length(rem_indices)==0)
    {
      break;
    }

    #Keep track of previous information
    degree_A <- main_degree_A[rem_indices];
    prev_jacc_df <- jacc_df;
    rm(jacc_df);
    gc();

  }
  return(ordered_list)
}

#Remove insignificant edges from graph and make the graph cleaner
prune_edges <- function(A,min_size_cluster)
{
  #Copy the sparse matrix
  B <- A;
  N <- nrow(A)
  k <- min_size_cluster
  for (i in 1:N)
  {
    flag <- 0;
    x_right <- 0;
    x_left <- 0;
    y_up <- 0;
    y_down <- 0;
    #For first k elements just check to right and to up
    if (i<=k)
    {
      x_right <- sum(A[i,(i:i+k-1)]);
      y_up <- sum(A[i:(i+k-1),i]);
      if (x_right<=0 && y_up<=0)
      {
        B[i,] <- 0;
        B[,i] <- 0;
      }
    }
    #For last k elements just check to left and to down
    else if (i>=N-k)
    {
      x_left <- sum(A[i,(i-k+1):i]);
      y_down <- sum(A[(i-k+1):i,i]);
      if (x_left<=0 && y_down<=0)
      {
        B[,i] <- 0;
        B[i,] <- 0;
      }
    }
    else
    {
      #Rest of the cases
      x_left <- sum(A[i,(i-k+1):i]);
      x_right <- sum(A[i,i:(i+k-1)]);
      y_up <- sum(A[i:(i+k-1),i]);
      y_down <- sum(A[(i-k+1):i,i])
      if (x_left>0 && x_right>0 && y_up>0 && y_down>0)
      {
        flag=1;
      }
      if (flag!=1)
      {
        B[,i] <- 0;
        B[i,] <- 0;
      }
    }
  }
  return(B);
}
