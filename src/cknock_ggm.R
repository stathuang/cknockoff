##################################################-
## Project: Conditional Knockoff
## Script purpose: Gaussian Graphical model
## Date: 2019-03-05
## Author: Dongming Huang (dhuang01@g.harvard.edu)
##################################################-

library(igraph)

# This function generates conditional knockoff for Gaussian graphical models.
# any one of the param. K, MaxK, scale 
# specify the number of folds used for data splitting
# param. threshold is the input n' for searching blocking sets
cknockoff.ggm = function(X, Adj, threshold, MaxK = Inf, K=NULL, 
                         scale = 2, method='mix') {
  n = nrow(X)
  p = ncol(X)
  
  # decide the number of data-splitting
  if(is.null(K)){
    K = floor(n / threshold / scale)
    if (K * scale <= 1) {
      print('blocking error(K==1)')
      return(NULL)
    }
    K = min(max(2, K), MaxK)
  }
  
  # the sample indices of the data-splitting
  sample.ind = matrix(c(sample(n), rep(0, ceiling(n / K) * K - n)), nr = K)
  
  Xk = c()  # initalize the knockoff
  blocked.count = rep(0, p)  # recording the times individual variables are blocked
  for (k in 1:K) {
    # sort by how often they are current blocked
    preempt.ind = order(blocked.count, decreasing = T)  
    # use the subroutine to generate conditional knockoff for this fold of data
    return2 = cknockoff.ggm.basic(
      X[sample.ind[k, ], ], 
      Adj = Adj, 
      var.ind = preempt.ind, 
      preempt.ind = NULL, 
      threshold = threshold, 
      method=method)
    
    blocked.count[return2$blocked.var] = blocked.count[return2$blocked.var] + 1
    Xk = rbind(Xk, return2$Xk)
  }
  
  if (sum(blocked.count==K) > 0) {
    print('Warning: trivial knockoff columns created.')
  }
  
  #reorder the sample indices
  permuted.ind = setdiff(c(t(sample.ind)), 0)
  Xk[order(permuted.ind), ]
}

# This function generates conditional knockoff for gaussian graphical model 
# with blocking set found by Graph_Expanding.ggm()
# Used as a subroutine for cknockoff.ggm()
# param. threshold is the n' input
# c(preempt.ind,var.ind ) is the $\pi$ input, the order to visit each variable
# Adj is the adjacent matrix of the graph
cknockoff.ggm.basic = function(X, Adj, var.ind = NULL, 
                               preempt.ind = NULL, threshold = NULL, 
                               method = 'sdp') {
  n = nrow(X)
  p = ncol(X)
  if (is.null(threshold))
    threshold = n
  if (threshold > n)
    threshold = n
  
  # search for a blocking set
  block = Graph_Expanding.ggm(
    Adj = Adj, 
    threshold = threshold, 
    var.ind = var.ind, 
    preempt.ind = preempt.ind
  )
  if (is.null(block)) {  # if no need to block any variable, use low-dimensional
    Xk = create.ldg(X)
  } else{
    # decompose the subgraph that deletes the blocking set
    graph = graph_from_adjacency_matrix(Adj[-block, -block] !=0, 
                                        mode = "undirected")
    membership = components(graph)$membership
    remaining = setdiff(1:p, block)
    
    Xk = X  # initialize the knockoff
    
    for (m in unique(membership)) {  # visit each connected component
      ind = remaining[membership == m]
      blocking.ind = block[which(colSums(Adj[ind, block, drop = F] != 0) > 0)]
      d = length(ind)
      b = length(blocking.ind) + 1
      if (n - b < 2 * d){
        print('Error: the graph is not n-separated by B')
        return(NULL)
      }
      
      # generate partial conditional knockoff for this component
      Xk[, ind] = cknofkoff.ldg.partial(X[, ind, drop = F], 
                                        X[, blocking.ind, drop = F], 
                                        method=method)
    }
  }
  return(list(Xk = Xk, blocked.var = block))
}


# greedy search for a block set by graph expanding
# as a subroutine for cknockoff.ggm.basic()
Graph_Expanding.ggm = function(Adj, threshold, var.ind = NULL, 
                               preempt.ind = NULL ) {
  p = ncol(Adj)
  if (is.null(var.ind)) {
    var.ind = 1:p
  }
  if (!is.null(preempt.ind)) {
    var.ind = c(preempt.ind, setdiff(var.ind, preempt.ind))
  }
  
  Adj = (Adj != 0)  # make sure the adjacent matrix is binary and symmetric
  if(sum(Adj!=t(Adj))){
    print('Warning: asymmetric adjacent matrix!')
  }
    
  blocked.var = c()  # initalize the blocking set
  for (i in 1:p) {
    ind = setdiff(which(Adj[var.ind[i] , ] != 0), var.ind[i])    # edge from i
    if (i > 1) {   # get the neighborhood of i in the expanded graph
      undrawn.neigh = setdiff(ind, var.ind[1:(i - 1)])
      regressor = c(ind, setdiff(intersect(ind, var.ind[1:(i - 1)]) , blocked.var))
    } else{
      undrawn.neigh = ind
      regressor = ind
    }
    if (length(regressor) + 3 > threshold) {   # if too many neighbors, block i
      blocked.var = c(blocked.var, var.ind[i])
    } else{
      Adj[undrawn.neigh, ind] = T  # undate the graph 
    }
  }
  return(blocked.var)
}


# [outdated]
# cknockoff.ggm.twofold = function(X, Adj, threshold, 
#     bAdaptive = F, method = 'sdp') {
#   # divide into 2 only
# 
# n = nrow(X)
# p = ncol(X)
# 
# halfsample = sample(n, floor(n / 2))
# Xk = X
# 
# 
# return1 = cknockoff.ggm.basic(X[halfsample, ], 
#     Adj = Adj, 
#     threshold = threshold, 
#     method = method)
# if (is.null(return1))
# return(NULL)
# return2 = cknockoff.ggm.basic(X[-halfsample, ], 
#     Adj = Adj, 
#     var.ind = 1:p, 
#     preempt.ind = return1$blocked.var, 
#     threshold = threshold, 
#     method = method)
# B0 = intersect(return1$blocked.var, return2$blocked.var)
# if (length(B0) > 0)
# print('Warning:trivial columns created.')
# 
# Xk[halfsample, ] = return1$Xk
# Xk[-halfsample, ] = return2$Xk
# Xk
# }
