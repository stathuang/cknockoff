# utility functions

# vector length
l2 = function(v)
  sqrt(sum(v ^ 2))

# Fast versions of diag(d) %*% X and X %*% diag(d).
`%diag*%` <- function(d, X)
  d * X  #diag(d) %*% X
`%*diag%` <- function(X, d)
  t(t(X) * d)


# scale the columns to have norm 1
normc = function(X, center = T) {
  X.centered = scale(X, center = center, scale = F)
  X.scaled = scale(X.centered, center = F, 
                   scale = sqrt(colSums(X.centered ^2)))
  X.scaled[, ] # No attributes
}

# selection errors, false positive and false negative
fp = function(selected, beta)
  sum(beta[selected] == 0) / max(1, length(selected))
fn = function (selected, beta)
  sum(beta[selected] != 0) / max(1, sum(beta != 0))


# create adjacent matrix for  spatial Graphical Models
SpatialGM = function(a, b) {
  p = a * b
  mat.adj = matrix(0, p, p)
  for (i.a in (1:a - 1)) {
    for (i.b in 1:b) {
      k = i.a * b + i.b
      if (i.a > 0)
        mat.adj[k, i.a * b - b + i.b] = 1
      if (i.a < a - 1)
        mat.adj[k, i.a * b + b + i.b] = 1
      if (i.b > 1)
        mat.adj[k, i.a * b + i.b - 1] = 1
      if (i.b < b)
        mat.adj[k, i.a * b + i.b + 1] = 1
    }
  }
  return(mat.adj)
}

Spatial3d = function(w, h, d) {
  p = w * h * d
  mat.adj = matrix(0, p, p)
  for (i.d in (1:d - 1)) {
    for (i.w in (1:w - 1)) {
      for (i.h in 1:h) {
        k = i.w * h + i.h + i.d * w * h
        if (i.w > 0)
          mat.adj[k, k - h] = 1
        if (i.w < w - 1)
          mat.adj[k, k + h] = 1
        if (i.h > 1)
          mat.adj[k, k - 1] = 1
        if (i.h < h)
          mat.adj[k, k + 1] = 1
        
        if (i.d > 0)
          mat.adj[k, k - w * h] = 1
        if (i.d < d - 1)
          mat.adj[k, k + w * h] = 1
      }
    }
  }
  return(mat.adj)
}


# used in sampling for ising model coupling from the past
sum.neighbors = function(v, S) {
  res = 0
  w = nrow(S)
  h = ncol(S)
  if (v[1] > 1)
    res = res + S[v[1] - 1, v[2]]
  if (v[1] < w)
    res = res + S[v[1] + 1, v[2]]
  if (v[2] > 1)
    res = res + S[v[1], v[2] - 1]
  if (v[2] < h)
    res = res + S[v[1], v[2] + 1]
}

# Exact Sampling for Ising model
CFPT.Ising = function(size = 50, temperature = 5) {
  beta = 1 / temperature
  
  Xmax = matrix(1, size, size)
  
  Xmin = -Xmax
  
  k = 1
  U = runif(1)
  V = matrix(sample(1:size, 2 ^ k , rep = T), nr = 2)
  
  while (sum(Xmax != Xmin) > 0 & k < 32) {
    V = cbind(matrix(sample(1:size, 2 ^ k , rep = T) , nr = 2), V)
    U = c(runif(2 ^ (k - 1)) ,  U)
    Xmax = matrix(1, size, size)
    
    Xmin = -Xmax
    
    for (i in 1:(2 ^ k)) {
      v = V[, i]
      u = U[i]
      Gmin = sum(sum.neighbors(v , Xmin))
      Gmax = sum(sum.neighbors(v, Xmax))
      
      if (u < 1 / (1 + exp(-beta * Gmin)))
      {
        Xmin[v[1], v[2]] = 1
      } else{
        Xmin[v[1], v[2]] = -1
      }
      
      if (u < 1 / (1 + exp(-beta * Gmax)))
      {
        Xmax[v[1], v[2]] = 1
      } else{
        Xmax[v[1], v[2]] = -1
      }
      
    }
    k = k + 1
  }
  
  if (k > 31)
    print('Not stop in 2^31 steps.')
  return(Xmax)
}
