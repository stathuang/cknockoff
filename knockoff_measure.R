## knockoff measures functions

mod.stat.glmnet_lambdadiff =  function (X, X_k, y, family = "gaussian", ...)
{
  res = stat.glmnet_lambdadiff(X, X_k, y, family = "gaussian", ...)
  zero.ind = colSums((X_k - X)^2)==0
  res[zero.ind]=0
  res
}

mod.stat.lasso_lambdadiff=function (X, X_k, y, ...) {
  res = stat.lasso_lambdadiff(X, X_k, y, ...)
  zero.ind = colSums((X_k - X)^2)==0
  res[zero.ind]=0
  res
}
mod.stat.lasso_lambdasmax =  function (X, X_k, y, ...)
{
  res = stat.lasso_lambdasmax(X, X_k, y, ...)
  zero.ind = colSums((X_k - X)^2)==0
  res[zero.ind]=0
  res
}
mod.stat.lasso_coefdiff =  function (X, X_k, y, ...)
{
  res = stat.lasso_coefdiff(X, X_k, y)#, ...)
  zero.ind = colSums((X_k - X)^2)< 1e-12 * length(y)
  res[zero.ind]=0
  res
}
mod.stat.stability_selection=function (X, X_k, y, ...)
{
  res = stat.stability_selection(X, X_k, y,...)
  zero.ind = colSums((X_k - X)^2)< 1e-12 * length(y)
  res[zero.ind]=0
  res
}

mod.stat.lasso_lambdadiff_bin=function (X, X_k, y, ...)
{
  res = stat.glmnet_lambdadiff(X, X_k, y, family = "binomial", ...)
  zero.ind = colSums((X_k - X)^2)==0
  res[zero.ind]=0
  res
}
mod.stat.lasso_lambdasmax_bin=function (X, X_k, y, ...)
{
  res = stat.lasso_lambdasmax_bin(X, X_k, y, family = "binomial", ...)
  zero.ind = colSums((X_k - X)^2)==0
  res[zero.ind]=0
  res
}
mod.stat.lasso_coefdiff_bin=function (X, X_k, y, cores = 2, ...)
{
  res=stat.lasso_coefdiff_bin(X, X_k, y,...)
  zero.ind = colSums((X_k - X)^2)< 1e-12 * length(y)
  res[zero.ind]=0
  res
}
mod.stat.random_forest=function (X, X_k, y, ...)
{
  res = stat.random_forest(X, X_k, y,...)
  zero.ind = colSums((X_k - X)^2)< 1e-12 * length(y)
  res[zero.ind]=0
  res
}


stat_scalreg = function(X, X_k, y,intercept=T) {
  if(intercept)
  {
    Z = scalreg(y = y, X = cbind(1, X, X_k))$coefficients[-1]
  }else{
    Z = scalreg(y = y, X = cbind(X, X_k))$coefficients
  }
  p = ncol(X)
  orig = 1:p
  res=abs(Z[orig]) - abs(Z[orig + p])
  zero.ind = colSums((X_k - X)^2)< 1e-12 * length(y)
  res[zero.ind]=0
  return(res)
}
