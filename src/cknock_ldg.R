##################################################-
## Project: Conditional Knockoff 
## Script purpose: Low-Dimensional Gaussian Models
## Date: 2019-03-05
## Author: Dongming Huang (dhuang01@g.harvard.edu)
##################################################-


# modify the decompose function in package `knockoff`
# return the SVD of a de-mean covariate matrix and compute the U matrix
decompose.demean=function (X)
{
  n = nrow(X)
  p = ncol(X)
  stopifnot(n >= 2 * p+1)
  result = knockoff:::canonical_svd(X)
  
  # differs from the fixed x knockoff
  Q = qr.Q(qr(cbind(1, result$u, knockoff:::rnorm_matrix(n, p))))
  u_perp = Q[, (p + 2):(2*p+1)] %*diag%sample(c(-1, 1), p, rep=T)
  
  result$u_perp = u_perp
  result
}

# Modify the create.solve_sdp() function in package `knockoff`.
# param. `scale.lower` is the constant in front of the s^{EQ}
# as the lower bound for s^{SDP}
create.solve_sdp_mix=function (Sigma, gaptol = 1e-06, maxit = 1000, 
                               scale.lower=0.1)
{
  stopifnot(isSymmetric(Sigma))
  G = cov2cor(Sigma)
  p = dim(G)[1]
  if (!knockoff:::is_posdef(G)) {
    stop("The covariance matrix is not positive-definite: cannot solve SDP", 
         immediate. = T)
  }
  
  # solve the equi-problem
  lambda_min = eigen(G, symmetric = T, only.values = T)$values[p]
  if (lambda_min < 0) {
    stop("In creation of equi-correlated knockoffs, while computing the smallest eigenvalue of the \n    covariance matrix. The covariance matrix is not positive-definite.")
  }
  s = rep(1, nrow(Sigma)) * min(2 * lambda_min, 1)
  psd = 0
  s_eps = 1e-08
  while (psd == 0) {
    psd = knockoff:::is_posdef(2 * G - diag(s * (1 - s_eps), length(s)))
    if (!psd) {
      s_eps = s_eps * 10
    }
  }
  s = s * (1 - s_eps) # this computes the s^{EQ}
  
  Cl1 = -rep(scale.lower*s[1], p)    # nagative the lower bound
  Al1 = -Matrix::Diagonal(p)
  Cl2 = rep(1, p)
  Al2 = Matrix::Diagonal(p)
  d_As = c(diag(p))
  As = Matrix::Diagonal(length(d_As), x = d_As)
  As = As[which(Matrix::rowSums(As) > 0), ]
  Cs = c(2 * G)
  A = cbind(Al1, Al2, As)
  C = matrix(c(Cl1, Cl2, Cs), 1)
  K = NULL
  K$s = p
  K$l = 2 * p
  b = rep(1, p)
  OPTIONS = NULL
  OPTIONS$gaptol = gaptol
  OPTIONS$maxit = maxit
  OPTIONS$logsummary = 0
  OPTIONS$outputstats = 0
  OPTIONS$print = 0
  sol = Rdsdp::dsdp(A, b, C, K, OPTIONS)
  if (!identical(sol$STATS$stype, "PDFeasible")) {
    warning("The SDP solver returned a non-feasible solution")
  }
  s = sol$y
  s[s < 0] = 0
  s[s > 1] = 1
  psd = 0
  s_eps = 1e-08
  while (psd == 0) {
    if (knockoff:::is_posdef(2 * G - diag(s * (1 - s_eps), length(s)))) {
      psd = 1
    }
    else {
      s_eps = s_eps * 10
    }
  }
  s = s * (1 - s_eps)
  if (max(s) == 0) {
    warning("In creation of SDP knockoffs, procedure failed. Knockoffs will have no power.", 
            immediate. = T)
  }
  return(s * diag(Sigma))
}

# Main function to generate conditional knocoff for n>2p
# param. `method` allows to choose different method to compute s
cknockoff.ldg=function (X, method='mix')
{
  n=nrow(X);p=ncol(X)
  X.center=X-matrix(colMeans(X), n, p, byrow=T) # de-mean
  X.svd = decompose.demean(X.center) # compute SVD and the U matrix
  tol = 1e-05
  d = X.svd$d
  d_inv = 1/d
  d_zeros = d <= tol * max(d)
  if (any(d_zeros)) {
    warning(paste("Data matrix is rank deficient.", "Model is not identifiable, but proceeding with SDP knockoffs"), 
            immediate. = T)
    d_inv[d_zeros] = 0
  }
  G = tcrossprod(X.svd$v %*diag% d)
  G_inv = tcrossprod(X.svd$v %*diag% d_inv)
  
  # different methods for computing s vector
  switch(method, 
         'sdp'={
           s = create.solve_sdp(G)
         }, 
         'equi'={
           s = create.solve_equi(G)
         }, 
         'asdp'={
           s = create.solve_asdp(G)
         }, 
         'mix'={
           s=create.solve_sdp_mix(G, scale.lower =0.1 )
         }, 
         'hybrid'={
           s=create.solve_sdp(G)*0.9+create.solve_equi(G)*0.1
         }
  )
  
  s[s <= tol] = 0
  C.svd = knockoff:::canonical_svd(2 * diag(s) - (s %diag*% G_inv %*diag%
                                                    s))
  X_ko = X - (X.center %*% G_inv %*diag% s) + (X.svd$u_perp %*diag%
                                                 sqrt(pmax(0, C.svd$d))) %*% t(C.svd$v)
}

# a subroutine used in cknofkoff.ldg.partial()
# for computing s and L based on the Gram matrix of R
create.ldg.block = function (R, U1, method = 'mix')
{
  p = ncol(R)
  G = crossprod(R)
  X.svd = eigen(G)
  
  tol = 1e-05
  d = sqrt(X.svd$values)
  d_inv = 1 / d
  d_zeros = d <= tol * max(d)
  if (any(d_zeros)) {
    warning(
      paste(
        "Data matrix is rank deficient.", 
        "Model is not identifiable, but proceeding with SDP knockoffs"
      ), 
      immediate. = T
    )
    d_inv[d_zeros] = 0
  }
  G_inv = tcrossprod(X.svd$vectors %*diag% d_inv)
  
  # different methods for computing s vector
  switch(
    method, 
    'sdp' = {
      s = create.solve_sdp(G)
    }, 
    'equi' = {
      s = create.solve_equi(G)
    }, 
    'asdp' = {
      s = create.solve_asdp(G)
    }, 
    'mix' = {
      s = create.solve_sdp_mix(G, scale.lower = 0.1)
    }, 
    'hybrid' = {
      s = create.solve_sdp(G) * 0.9 + create.solve_equi(G) * 0.1
    }
  )
  
  s[s <= tol] = 0
  C.svd = knockoff:::canonical_svd(2 * diag(x = s) - (s %diag*% G_inv %*diag%
                                                        s))
  X_ko = R - (R %*% G_inv %*diag% s) + (U1 %*diag%
                                          sqrt(pmax(0, C.svd$d))) %*% t(C.svd$v)
}

# Generating Parital Knockoff for Gaussian models
cknofkoff.ldg.partial=function(Xv, Xb, method='mix'){
  d=ncol(Xv)
  b=ncol(Xb)+1
  n=nrow(Xv)
  if(nrow(Xb)!=n){
    print('Error: row numbers not match!')
  }
  #Gram-Schmidt Orthonormalization; the sign need adjustment
  GS = qr.Q(qr((cbind(
    1, Xb, Xv, matrix(rnorm(n * d), nr =n)
  ))))
  U1 = GS[, -(1:(b + d)), drop = F] %*diag% sample(c(-1, 1), d, rep = T)
  
  proj.mat = GS[, 1:b]  # projection matrix on to the columns of Xb
  projected = tcrossprod(proj.mat) %*% Xv  # projection of Xv 
  R = Xv - projected  # residuals
  if (ncol(R) == 1) {
    Xkv = projected + l2(R) * U1  # if V is univariate, rescaling U1 is enough
  } else{
    Xkv = projected + create.ldg.block(R, U1, method)  # compute s and L
  }
  return(Xkv)
}

