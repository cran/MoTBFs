# library(bnlearn)
# moral
# bnlearn:::dag2ug.backend
# 
# 
# dag = model2network("[A][C][F][B|A][D|A:C][E|B:F]")
# 
# graphviz.plot(dag)
# mdag = moral(dag)
# graphviz.plot(mdag)
# 
# parents(mdag, "B")
# children(mdag, "B")
# 
# nbr(mdag, "B")
# 
# elimin_order = c()
# graphviz.plot(mdag)
# len = sapply(lapply(nodes(mdag), nbr, x = mdag), length)
# names(len) = nodes(mdag)
# len
# which(len == min(len))
# rm.node = names(which.min(len))
# elimin_order = c(elimin_order, rm.node)
# mdag = remove.node(mdag, rm.node)
# 
# noHidden = c("G","C")
# min_nbr = function(dag, noHidden = NULL){
#   # moralize and triangulate DAG
#   mdag = moral(dag)
#   
#   # initialize elimination ordering
#   elimin_order = c()
#   
#   # vector of hidden node names
#   n = nodes(mdag)
#   n = n[!(n %in%noHidden)]
#   
#   while(length(n)>1){
#     # compute number of neighbors of each node
#     len = sapply(lapply(n, nbr, x = mdag), length)
#     names(len) = n
#     
#     # find node with min nbr
#     rm.node = names(which.min(len))
#     
#     # store node
#     elimin_order = c(elimin_order, rm.node)
#     
#     # remove node from moralized DAG
#     # graphviz.plot(mdag)
#     mdag = remove.node(mdag, rm.node)
#     
#     # update vector of nodes
#     n = nodes(mdag)
#     n = n[!(n %in%noHidden)]
#   }
#   
#   # add last one node to the ordering vector
#   elimin_order = c(elimin_order, n)
#   return(elimin_order)
# }
# 
# 
# dag = model2network("[A][B][E][G][C|A:B][D|B][F|A:D:E:G]")
# graphviz.plot(dag)
# graphviz.plot(moral(dag))
# min_nbr(dag)
# 
# 
# 
# dag = model2network("[A][B|A][C|A][D|B:C][E|C]")
# graphviz.plot(dag)
# mdag = moral(dag)
# dag2 = remove.node(dag, "B")
# graphviz.plot(mdag)
# 
# degree(dag, "C")
# m = undirected.arcs(mdag)
# mdag$arcs = rbind(m, c('A','E'))
# 
# 
# 
# 
# 
# 
# 
# 
# noHidden = NULL

# Function to implement the Min-Degree Heuristic
#' @importFrom bnlearn moral nodes amat
#' @noRd
min_degree_heuristic <- function(dag, noHidden = NULL) {
  
  adj_matrix = amat(moral(dag))
  
  # Number of variables (nodes)
  num_vars <- nrow(adj_matrix)
  
  # List to store the elimination order
  elimination_order <- c()
  
  # Keep track of remaining variables
  remaining_vars = 1:num_vars
  # except for target and evidence variables
  remaining_vars = setdiff(remaining_vars,which(nodes(dag)%in%noHidden))
  
  # Repeat until all variables are eliminated
  while (length(remaining_vars) > 1) {
    
    # Compute degrees of all remaining variables
    degrees <- colSums(adj_matrix[remaining_vars, remaining_vars])
    
    # Select the variable with the smallest degree
    min_degree_var <- remaining_vars[which.min(degrees)]
    
    # Add the variable to the elimination order
    elimination_order <- c(elimination_order, min_degree_var)
    
    # Get the neighbors of the selected variable
    neighbors <- which(adj_matrix[min_degree_var, ] == 1)
    neighbors <- neighbors[neighbors %in% remaining_vars]  # Only keep remaining variables
    
    # Add fill-in edges between all pairs of neighbors
    if (length(neighbors) > 1) {
      for (i in 1:(length(neighbors) - 1)) {
        for (j in (i + 1):length(neighbors)) {
          adj_matrix[neighbors[i], neighbors[j]] <- 1
          adj_matrix[neighbors[j], neighbors[i]] <- 1
        }
      }
    }
    
    # Eliminate the variable by removing it from the adjacency matrix
    remaining_vars <- setdiff(remaining_vars, min_degree_var)
    adj_matrix[min_degree_var, ] <- 0
    adj_matrix[, min_degree_var] <- 0
  }
  elimination_order = colnames(adj_matrix)[c(elimination_order,remaining_vars)]
  # Return the elimination order
  return(elimination_order)
}








# Function to implement the Min-Fill Heuristic
#' @importFrom bnlearn moral nodes amat
#' @noRd
min_fill_heuristic <- function(dag, noHidden = NULL) {
  
  adj_matrix = amat(moral(dag))
  
  # Number of variables (nodes)
  num_vars <- nrow(adj_matrix)
  
  # List to store the elimination order
  elimination_order <- c()
  
  # Keep track of remaining variables
  remaining_vars <- 1:num_vars
  
  # except for target and evidence variables
  remaining_vars = setdiff(remaining_vars,which(nodes(dag)%in%noHidden))
  
  # Repeat until all variables are eliminated
  while (length(remaining_vars) > 1) {
    
    # Initialize a variable to store the number of fill-ins for each variable
    fill_in_count <- rep(0, length(remaining_vars))
    
    # For each remaining variable, calculate the number of fill-ins
    for (i in seq_along(remaining_vars)) {
      var <- remaining_vars[i]
      
      # Get the neighbors of the current variable
      neighbors <- which(adj_matrix[var, ] == 1)
      neighbors <- neighbors[neighbors %in% remaining_vars]  # Only consider remaining variables
      
      # Count how many pairs of neighbors are not connected (fill-ins needed)
      if (length(neighbors) > 1) {
        for (j in 1:(length(neighbors) - 1)) {
          for (k in (j + 1):length(neighbors)) {
            if (adj_matrix[neighbors[j], neighbors[k]] == 0) {
              fill_in_count[i] <- fill_in_count[i] + 1
            }
          }
        }
      }
    }
    
    # Select the variable with the smallest number of fill-ins
    min_fill_var <- remaining_vars[which.min(fill_in_count)]
    
    # Add the variable to the elimination order
    elimination_order <- c(elimination_order, min_fill_var)
    
    # Get the neighbors of the selected variable
    neighbors <- which(adj_matrix[min_fill_var, ] == 1)
    neighbors <- neighbors[neighbors %in% remaining_vars]  # Only keep remaining variables
    
    # Add fill-in edges between all pairs of neighbors
    if (length(neighbors) > 1) {
      for (i in 1:(length(neighbors) - 1)) {
        for (j in (i + 1):length(neighbors)) {
          adj_matrix[neighbors[i], neighbors[j]] <- 1
          adj_matrix[neighbors[j], neighbors[i]] <- 1
        }
      }
    }
    
    # Eliminate the variable by removing it from the adjacency matrix
    remaining_vars <- setdiff(remaining_vars, min_fill_var)
    adj_matrix[min_fill_var, ] <- 0
    adj_matrix[, min_fill_var] <- 0
  }
  elimination_order = colnames(adj_matrix)[c(elimination_order,remaining_vars)]
  # Return the elimination order
  return(elimination_order)
}

