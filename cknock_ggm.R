## gaussian graphical model

library(igraph)

# search for a block set
# Alg3
Graph_Expanding.ggm = function(Adj,  threshold,  var.ind = NULL,
                             preempt.ind = NULL ) {
  p = ncol(Adj)
  if (is.null(var.ind)) {
    var.ind = 1:p
  }

  if (!is.null(preempt.ind)) {
    var.ind = c(preempt.ind, setdiff(var.ind, preempt.ind))
  }

  Adj = (Adj != 0)

  blocked.var = c()
  for (i in 1:p) {
    ind = setdiff(which(Adj[var.ind[i] , ] != 0), var.ind[i])         # out edge
    if (i > 1) {
      undrawn.neigh = setdiff(ind, var.ind[1:(i - 1)])
      regressor = c(ind, setdiff(intersect(ind, var.ind[1:(i - 1)]) , blocked.var))
    } else{
      undrawn.neigh = ind
      regressor = ind
    }
    if (length(regressor) + 3 > threshold) { # if too much neighbors, block it
      blocked.var = c(blocked.var, var.ind[i])
    } else{
      Adj[undrawn.neigh, ind] = T
    }
  }
  return(blocked.var)
}



cknockoff.ggm.basic = function(X, Adj, var.ind = NULL,
                                          preempt.ind = NULL, threshold = NULL,
                                          method = 'sdp') {
    n = nrow(X)
    p = ncol(X)
    if (is.null(threshold))
      threshold = n
    if (threshold > n)
      threshold = n

    block = Graph_Expanding.ggm(
      Adj = Adj,
      threshold = threshold,
      var.ind = var.ind,
      preempt.ind = preempt.ind
    )
    if (is.null(block)) {
      Xk = create.ldg(X)
    } else{
      graph = graph_from_adjacency_matrix(Adj[-block, -block] !=0,
                                          mode = "undirected")
      membership = components(graph)$membership
      remaining = setdiff(1:p, block)
      Xk = X

      for (m in unique(membership)) {
        ind = remaining[membership == m]

        blocking.ind = block[which(colSums(Adj[ind, block, drop = F] != 0) > 0)]
        d = length(ind)
        b = length(blocking.ind) + 1
        if (n - b < 2 * d){
          print('Error: the graph is not n-separated by B')
          return(NULL)
        }

          Xk[, ind] = cknofkoff.ldg.partial(X[, blocking.ind, drop = F], X[, ind, drop = F],method=method)
      }
    }
    return(list(Xk = Xk, blocked.var = block))
  }


cknockoff.ggm = function(X, Adj, threshold, MaxK = Inf, scale = 2) {
  n = nrow(X)
  p = ncol(X)
  K = floor(n / threshold / scale)

  if (K * scale <= 1) {
    print('blocking error(K==1)')
    return(NULL)
  }

  K = min(max(2, K), MaxK)
  sample.ind = matrix(c(sample(n), rep(0, ceiling(n / K) * K - n)), nr = K)

  Xk = c()
  blocked.count = rep(0, p)
  blocking = rep(T, p)
  for (k in 1:K) {
    preempt.ind = order(blocked.count, decreasing = T)
    return2 = cknockoff.ggm.basic(
      X[sample.ind[k,],],
      Adj = Adj,
      var.ind = preempt.ind,
      preempt.ind = NULL,
      threshold = threshold
    )
    blocking[setdiff(1:p, return2$blocked.var)] = F
    blocked.count[return2$blocked.var] = blocked.count[return2$blocked.var] + 1
    Xk = rbind(Xk, return2$Xk)
  }

  if (sum(blocking) > 0) {
    print('Warning:trivial columns created.')
  }

  #reorder the sample id!!
  permuted.ind = setdiff(c(t(sample.ind)), 0)
  Xk[order(permuted.ind),]
}

cknockoff.ggm.twofold = function(X, Adj, threshold, bAdaptive = F, method = 'sdp') {
  # divide into  2 only

  n = nrow(X)
  p = ncol(X)

  halfsample = sample(n, floor(n / 2))
  Xk = X


  return1 = cknockoff.ggm.basic(X[halfsample, ],
                                Adj = Adj,
                                threshold = threshold,
                                method = method)
  if (is.null(return1))
    return(NULL)
  return2 = cknockoff.ggm.basic(X[-halfsample, ],
                                Adj = Adj,
                                var.ind = 1:p,
                                preempt.ind = return1$blocked.var,
                                threshold = threshold,
                                method = method)
  B0 = intersect(return1$blocked.var, return2$blocked.var)
  if (length(B0) > 0)
    print('Warning:trivial columns created.')

  Xk[halfsample,] = return1$Xk
  Xk[-halfsample,] = return2$Xk
  Xk
}
