
# Conditional Knockoff for Low-Dimensional Gaussian Models

# modify the decompose function in package `knockoof`
decompose.demean=function (X)
{
  n = nrow(X)
  p = ncol(X)
  stopifnot(n >= 2 * p+1)
  result = knockoff:::canonical_svd(X)

  # differs from the fixed x knockoff
  Q = qr.Q(qr(cbind(1, result$u, knockoff:::rnorm_matrix(n,p))))
  u_perp = Q[, (p + 2):(2*p+1)]  %*diag%sample(c(-1,1),p,rep=T)

  result$u_perp = u_perp
  result
}

## modify the create.solve_sdp function in package `knockoff`
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
    stop("In creation of equi-correlated knockoffs, while computing the smallest eigenvalue of the \n                covariance matrix. The covariance matrix is not positive-definite.")
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
  s = s * (1 - s_eps)

  Cl1 = -rep(scale.lower*s[1], p) # nagative the lower bound
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

cknockoff.ldg=function (X,method='mix')
{
  n=nrow(X);p=ncol(X)
  X.center=X-matrix(colMeans(X),n,p,byrow=T)

  X.svd = decompose.demean(X.center)
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
           s=create.solve_sdp_mix(G,scale.lower =0.1 )
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

## partial conditional knockoff for ldg

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

cknofkoff.ldg.partial=function(Xv,Xb,method='mix'){
  d=ncol(Xv)
  b=ncol(Xb)+1
  n=nrow(Xv)
  if(nrow(Xb)!=n){
    print('Error: row numbers not match!')
  }
  #Gram-Schmidt Orthonormalization; the sign need adjustment
  GS = qr.Q(qr((cbind(
    1, Xv, Xb,  matrix(rnorm(n * d), nr =n)
  ))))
  U1 = GS[, -(1:(b + d)), drop = F] %*diag% sample(c(-1, 1), d, rep = T)
  proj.mat = GS[, 1:b]
  projected = tcrossprod(proj.mat) %*% Xv
  R = Xv - projected
  if (ncol(R) == 1) {
    Xkv = projected + l2(R) * U1
  } else{
    Xkv = projected + create.ldg.block(R, U1, method)
  }
  return(Xkv)
}

