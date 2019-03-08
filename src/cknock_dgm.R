##################################################-
## Project: Conditional Knockoff
## Script purpose: General Discrete Graphical Models
## Date: 2019-03-05
## Author: Dongming Huang (dhuang01@g.harvard.edu)
## Status: Finished commenting up to line 103
##################################################-

# General conditional knockoff for discrete graphical models----

# a subroutine used by cknofkoff.DG.blocking.basic()
# randomly permute the entries of input `target` in each categories 
# defined by the combination of value of the input `fixed`. 
Restricted_Knockoff=function(target, fixed,K=2){
  n=length(target)
  if( nrow(fixed)!=n ){
    print('dimension not match!')
    return(NULL)
  }

  if(length(K)==1){K=rep(K,ncol(fixed) )} # the number of states for the fixed variables
  
  # divide the data by the combination of `fixed`
  group.ind=split(1:n,  lapply(1:ncol(fixed),function(v){
    factor(fixed[,v],levels=0:(K[v]-1))
  }))
  lens=unlist(lapply(group.ind,length))
  
  xk=rep(0,n)  # initialize the knockoff variable
  for( ii in 1:length(group.ind)){
    if(lens[ii]==0)next

    ind=group.ind[[ii]]
    if(lens[ii]==1){
      xk[ ind] = target[ind]  # if this category only has one unit
    }else{
      xk[ind] = sample(target[ind])  # randomly permute
    }
  }
  xk
}

# a subroutine used in cknofkoff.DG.blocking()
# generate conditional knockoff for variables in `var.ind` in 
# discrete graphical model with a cut set, i.e. the complement of `var.ind`
cknofkoff.DG.blocking.basic = function(x, Graph, var.ind,K=2) {
  for( j in var.ind){
    neighbour=setdiff(which(Graph[j,]!=0),j)
    if( length(intersect(neighbour,var.ind)  ) >0){
      print(j)
      print( ' The input variables are connected!')
      return(NULL)
    }
  }

  n=nrow(x)
  p=ncol(x)
  xk=x
  for( j in var.ind){
    neighbour=setdiff(which(Graph[j,]!=0),j)
    xk[,j]= Restricted_Knockoff(x[,j],x[,neighbour,drop=F],K=K)

  }
  return(xk)
}

# Generate conditional knockoff for discrete graphical model.
# `Graph`` is the adjacent matrix.
# `comps.ind` is a list of variables indices, each of denote a cut set.
# if only split data into 2 folds, `comps.ind` can include only one cut set.
# K is the number of states of each column in `fixed`, 
# which can be a scaler if they are equal. 
# individual variable should take value from 0 to K-1.
cknofkoff.DG.blocking=function(x,Graph,comps.ind,K=2){ 
  x=as.matrix(x)  # avoid using data.frame class
  
  bMinus=F  # check if the data are binary and coded as +1 and -1
  if(min(x)<0 & K==2){
    x[x==-1]=0
    bMinus=T
  }
  n=nrow(x);p=ncol(x)
  
  xk=x  # initialize the knockoff

  # if input a vector as a cut set, then the rest makes up the other cut set 
  if(is.vector(comps.ind)){ 
    n.comp=2
    comps.ind=list(comps.ind,setdiff(1:p,comps.ind))
  }

  # data splitting
  mat.Ind=matrix(c(sample(n),rep(0,ceiling(n/n.comp)*n.comp-n)), nr=n.comp)
  for( i.comp in 1:length(comps.ind) ){
    sample.ind= mat.Ind[i.comp,]
    # generate conditional knockoff for the i-th fold of data
    xk[sample.ind,]=cknofkoff.DG.blocking.basic(x[sample.ind,,drop=F],  
                                      Graph,comps.ind[[i.comp]],K=K)
  }

  if(bMinus)xk[xk==0]=-1  # change back to the original coding
  return(xk)
}

# Graph-expanding-------
cknofkoff.DG.blocking.basic.enhanced = function(x, Graph,var.ind,K=2) {
  expandedgraph=Graph
  for( j in var.ind){
    neighbour=setdiff(which(Graph[j,]!=0),j)
    expandedgraph[neighbour,neighbour]=1
    expandedgraph=rbind(expandedgraph,expandedgraph[j,])
    expandedgraph=cbind(expandedgraph,c(expandedgraph[j,],1))
    if( length(intersect(neighbour,var.ind)  ) >0){
      print(j)
      print( ' The input variables are connected!')
      return(NULL)
    }
  }

  n=nrow(x)
  p=ncol(x)
  xk=x
  for( j in var.ind){
    neighbour=setdiff(which(Graph[j,]!=0),j)
    xk[,j]= Restricted_Knockoff(x[,j],x[,neighbour,drop=F],K=K)

  }
  return(list(xk=xk,newgraph=expandedgraph))
}


cknofkoff.DG.Enhanced.blocking=function(x,Graph,comps.ind,K=2){
  x=as.matrix(x)
  bMinus=F
  if(min(x)<0 & K==2){
    x[x==-1]=0
    bMinus=T
  }


  n=nrow(x);p=ncol(x)


  n.comp=8

  mat.Ind=matrix(c(sample(n),rep(0,ceiling(n/n.comp)*n.comp-n)), nr=n.comp)
  xk=x

  for( i.comp in 1:8 ){

    sample.ind= mat.Ind[i.comp,]

    var.ind=unlist(comps.ind[i.comp%%2+1+2*(0:3)])

    input.x=x[sample.ind,,drop=F]
    output1=cknofkoff.DG.blocking.basic.enhanced(input.x,
                                       Graph,var.ind,K=2)
    output.xk=output1$xk

    B=setdiff(1:p,var.ind)
    var.ind2=comps.ind[[i.comp]]
    output.xk2=cknofkoff.DG.blocking.basic(cbind(input.x,output.xk[,var.ind,drop=F]  ),
                                 output1$newgraph,var.ind2   , K=2)
    xk[sample.ind,]=output.xk
    xk[sample.ind,var.ind2]=output.xk2[,var.ind2,drop=F]

  }
  if(bMinus)xk[xk==0]=-1
  return(xk)
}


# Improved Blocking for Markov chain--------
# Generate a exchangable contingency table used in  Alg 15
Gibbs.ContingencyTable=function(S,len=10){
  Sk=S
  dims=dim(S)
  for(l in 1:len)  {
    proposal=Sk

    Rs=sample(dims[1],2,rep=F)
    Cs=sample(dims[2],2,rep=F)
    Ds=sample(dims[3],2,rep=F)

    for(ir in 1:2){
      for(ic in 1:2){
        for(id in 1:2){
          proposal[Rs[ir],Cs[ic],Ds[id]]=proposal[Rs[ir],Cs[ic],Ds[id]]+(-1)^(ir+ic+id)
        }
      }
    }

    if(min(proposal)<0){
      next
    }
    alpha=1
    for(ir in 1:2){
      for(ic in 1:2){
        if(  (ir+ic) %%2 ==0) {
          alpha=alpha*Sk[Rs[ir],Cs[ic],Ds[1]]/( Sk[Rs[ir],Cs[ic],Ds[2]] +1 )
        }else{
          alpha=alpha*Sk[Rs[ir],Cs[ic],Ds[2]]/( Sk[Rs[ir],Cs[ic],Ds[1]] +1 )
        }
      }
    }
    if(runif(1)<alpha)Sk=proposal
  }
  return(Sk)
}


cknofkoff.MC.Restricted_Knockoff2=function(xa,xb,target,K=2,n.MH=10){# this algorithm only works for binary
  n=length(target)
  if( length(xb)!=n |length(xa)!=n){
    return(NULL)
  }

  S=table(factor(xa,levels=0:(K-1)),
          factor(xb,levels=0:(K-1)),
          factor(target,levels=0:(K-1)) )
  # set.seed(10)
  Sk=Gibbs.ContingencyTable(S,len=n.MH)

  xk=rep(0,n)
  taba=as.numeric(as.character(dimnames(S)[[1]]))
  tabb=as.numeric(as.character(dimnames(S)[[2]]))
  tabc=as.numeric(as.character(dimnames(S)[[3]]))

  for( ir in seq_len(length(taba)))
    for(ic in seq_len(length(tabb))){
      ind=which(xa==taba[ir]&xb==tabb[ic])
      if(length(ind)==0)next
      if(length(ind)==1){
        xk[ ind] =rep(tabc,times=Sk[ir,ic,])
      }else{
        xk[ ind] = sample(rep(tabc,times=Sk[ir,ic,]))
      }
    }
  xk
}

cknofkoff.MC.Restricted_Knockoff1=function(xa,target,K=2){# this algorithm only works for binary
  n=length(target)
  if( length(xa)!=n){
    return(NULL)
  }


  xk=rep(0,n)
  taba=0:(K-1)

  for( ir in seq_len(length(taba))){
    ind=which(xa==taba[ir])
    if(length(ind)==0)next
    if(length(ind)==1){
      xk[ ind] = target[ind]
    }else{
      xk[ ind] = sample(target[ind])
    }
  }
  xk
}


cknofkoff.MC.blocking.basic = function(x,var.ind,K=2,order=1,n.MH=10) {

  # check the sub-chain is separaterd
  if(min(diff( sort(var.ind,decreasing = F) )) < order){
    print('the input sub-chain is not separated!')
    return(NULL)
  }


  n=nrow(x)
  p=ncol(x)
  xk=x
  for( j in var.ind){

    if(j>1 & j<p){
      xk[,j]=cknofkoff.MC.Restricted_Knockoff2(x[,j-1],x[,j+1],x[,j],K=K,n.MH=n.MH)
    }else if(j==1){
      xk[,j]=cknofkoff.MC.Restricted_Knockoff1(x[,j+1],x[,j],K=K)
    }else {
      xk[,j]=cknofkoff.MC.Restricted_Knockoff1(x[,j-1],x[,j],K=K)
    }

  }
  return(xk)
}


cknofkoff.MC.blocking=function(x,K=2,n.MH=10){
  x=as.matrix(x)  # this single line double the speed!!! # data.frame is much slower than matrx
  bMinus=F
  if(min(x)<0 & K==2){
    x[x==-1]=0
    bMinus=T
  }

  n=nrow(x);p=ncol(x)
  sample.ind=sample(n,floor(n/2))
  cov.ind=(1:floor((1+p)/2))*2-1
  xk=x
  xk[sample.ind,]=cknofkoff.MC.blocking.basic(x[sample.ind,,drop=F],cov.ind,K=K,n.MH=n.MH)
  xk[-(sample.ind),]=cknofkoff.MC.blocking.basic(x[-sample.ind,,drop=F],setdiff(1:p,cov.ind),K=K,n.MH=n.MH)
  if(bMinus)xk[xk==0]=-1
  return(xk)
}

cknofkoff.MC.Enhanced.blocking=function(x,Graph,comps.ind,nFold=2,K=2){
  x=as.matrix(x)
  bMinus=F
  if(min(x)<0 & K==2){
    x[x==-1]=0
    bMinus=T
  }

  n=nrow(x);p=ncol(x)

  if(nFold==2){
    n.comp=2
    comps.ind=list(comps.ind,setdiff(1:p,comps.ind))

    mat.Ind=matrix(c(sample(n),rep(0,ceiling(n/n.comp)*n.comp-n)), nr=n.comp)
    xk=x

    for( i.comp in 1:nFold){
      sample.ind= mat.Ind[i.comp,]
      var.ind=comps.ind[[i.comp]]
      input.x=x[sample.ind,,drop=F]
      output1=cknofkoff.DG.blocking.basic.enhanced(input.x,
                                         Graph,var.ind,K=2)
      output.xk=output1$xk


      B=setdiff(1:p,var.ind)
      var.ind2=B[2*(1:floor(length(B)/2))]
      output.xk2=cknofkoff.DG.blocking.basic(cbind(input.x,output.xk[,var.ind,drop=F]  ),
                                   output1$newgraph,var.ind2   , K=2)
      xk[sample.ind,]=output.xk
      xk[sample.ind,var.ind2]=output.xk2[,var.ind2,drop=F]
    }
  }
  else if(nFold==4){
    n.comp=4
    comps.ind=list(comps.ind,setdiff(1:p,comps.ind))

    mat.Ind=matrix(c(sample(n),rep(0,ceiling(n/n.comp)*n.comp-n)), nr=n.comp)
    xk=x
    for( i.fold in 1:nFold ){
      i.comp=(i.fold-1)%%2+1
      r.comp=(i.fold-1)%/%2
      var.ind=comps.ind[[i.comp]]


      sample.ind= mat.Ind[i.fold,]
      input.x=x[sample.ind,,drop=F]

      output1=cknofkoff.DG.blocking.basic.enhanced(input.x,Graph,var.ind,K=2)
      output.xk=output1$xk

      B=setdiff(1:p,var.ind)
      var.ind2=B[2*(1:floor((r.comp+length(B))/2)) - r.comp]
      output.xk2=cknofkoff.DG.blocking.basic(cbind(input.x,output.xk[,var.ind,drop=F]  ),
                                   output1$newgraph,var.ind2   , K=2)
      xk[sample.ind,]=output.xk
      xk[sample.ind,var.ind2]=output.xk2[,var.ind2,drop=F]


    }


  }

  if(bMinus)xk[xk==0]=-1
  return(xk)
}


### SCIP--------

# Checking the pairwise contingency table
CheckMatch = function(x1, x2, xk1, xk2,K=2) {
  for( i in 1:(K-1))
    for( j in 1:(K-1))
    {
      if(sum((x1==i) & (x2==j)) != sum((xk1==i) & (xk2==j) ))
        return(F)
    }
  return(T)
}
# Checking the marginal frequency
# Q:takes up 47% running time?
CheckMarginMatch=function(x1,xk1,K=2){
  for( l in 1:(K-1)){
    if(sum(x1==l) !=sum(xk1==l))
      return(F)
  }
  return (T)
}

BinToInt = function(v,n,K=2) {
  ans = sum(K ^ (1:n - 1) * v[1:n])
  if (ans == 0)
    ans = K ^ n
  ans
}

ItoB2=function(v,n)as.numeric(intToBits(v))[1:n]

ItoB=function(v,n,K=2){
  if(K==2)return(as.numeric(intToBits(v))[1:n])
  res=v
  conv=rep(0,n)
  for(i in 1:n){
    conv[i] = res %% K
    res= res %/% K
    if(res==0)break
  }
  conv
}
cknofkoff.MC.SCIP = function(x,J,K) {
  x=as.matrix(x)  # this single line doubles the speed!!! # data.frame is much slower than matrx
  bMinus=F
  if(min(x)<0 & K==2){
    x[x==-1]=0
    bMinus=T
  }

  n=nrow(x)
  N=K^n
  qx = apply(x, 2, BinToInt,n=n,K=K)
  # create knockoff
  Trans = array(dim = c(J, K ^ n, K ^ n), data = 0)
  xk = matrix(0, n, J)
  qxk = rep(0, J) # the compact representation of xk

  for (j in 1:J) {
    temp.Trans = array(dim = c(K ^ n, K ^ n), data = 0)
    NonZeroCols=rep(F, K^n)

    if (j == 1) {
      for (qj in 1:N) {
        qj.x = ItoB(qj,n,K)
        if( !CheckMarginMatch(qj.x,x[, j],K))
          next

        if (j < J) {
          for (qk in 1:N) {
            qk.x = ItoB(qk,n,K)
            if( !CheckMarginMatch(qk.x,x[, j+1],K) )
              next
            if (!CheckMatch(x[, j], x[, j + 1], qj.x,  qk.x,K))
              next
            temp.Trans[qj, qk] = 1
            NonZeroCols[qk]=T
          }
        } else{
          temp.Trans[qj, 1] = 1
          NonZeroCols[1]=T
        }

      }

      for( qk in which(NonZeroCols)){
        Trans[1,,qk]=temp.Trans[,qk] / sum(temp.Trans[,qk])
      }

    } else{
      qi = qxk[j - 1]
      qi.x = xk[, j - 1]

      for (qj in 1:N) {
        qj.x = ItoB(qj,n,K)

        if( !CheckMarginMatch(qj.x,x[, j],K) )
          next
        if (!CheckMatch(x[, j], x[, j - 1], qj.x,  qi.x,K) |
            !CheckMatch(x[, j], x[, j - 1], qj.x,  x[, j - 1],K))
          next

        if (j < J) {
          for (qk in 1:N) {
            qk.x = ItoB(qk,n,K)

            if( !CheckMarginMatch(qk.x,x[, j + 1],K) )
              next
            if (!CheckMatch(x[, j], x[, j + 1], qj.x,  qk.x,K))
              next
            temp.Trans[qj, qk] = 1
            NonZeroCols[qk]=T
          }
        } else{
          temp.Trans[qj, 1] = 1
          NonZeroCols[1]=T
        }

      }

      for( qk in which(NonZeroCols)){
        v=temp.Trans[,qk] * Trans[j - 1, qi ,]
        Trans[j,,qk]=v/sum(v)
      }
    }

    # draw a sample
    if (j == J) {
      qxk[j] = sample(N, 1, prob = Trans[j, , 1])
    } else{
      qxk[j] = sample(N, 1, prob = Trans[j, , qx[j + 1]])
    }
    xk[, j] = ItoB(qxk[j],n,K)

  }

  if(bMinus)xk[xk==0]=-1
  return(xk)
}


cknofkoff.MC.SCIP.sparse = function(x,J,K) { # the function itself takes 69% time
  x=as.matrix(x)
  bMinus=F
  if(min(x)<0 & K==2){
    x[x==-1]=0
    bMinus=T
  }

  n=nrow(x)
  N=K^n
  qx = apply(x, 2, BinToInt,n=n,K=K)
  # create knockoff
  lastTrans = matrix(data=0, K ^ n, K ^ n) #Matrix(data=0, K ^ n, K ^ n, sparse=T)
  xk = matrix(0, n, J)
  qxk = rep(0, J) # the compact representation of xk

  for (j in 1:J) {
    temp.Trans =matrix(data=0, K ^ n, K ^ n) #  Matrix(data=0, K ^ n, K ^ n, sparse=T)
    NonZeroCols=rep(F, K^n)

    if (j == 1) {
      for (qj in 1:N) {
        qj.x = ItoB(qj,n,K)
        # if (sum(qj.x) != sum(x[, j]))
        #   next
        if( !CheckMarginMatch(qj.x,x[, j],K)) # all CheckMarginMatch() only take 1.6% time
          next

        if (j < J) {
          for (qk in 1:N) {
            qk.x = ItoB(qk,n,K)
            # if (sum(qk.x) != sum(x[, j + 1]))
            #   next
            if( !CheckMarginMatch(qk.x,x[, j+1],K) )
              next
            if (!CheckMatch(x[, j], x[, j + 1], qj.x,  qk.x,K))
              next
            Debug.Inner(n)
            temp.Trans[qj, qk] = 1
            NonZeroCols[qk]=T
          }
        } else{
          Debug.Inner(n)
          # for j==J, no qk
          temp.Trans[qj, 1] = 1
          NonZeroCols[1]=T
        }

      }

      # inefficient
      # temp.Trans <-
      #   apply(temp.Trans  , 2, function(v)    # inefficient!!
      #     if(sum(v) > 0)
      #       v / sum(v)
      #     else
      #       0 * v)  # apply(,2,)keep the direction of the matrix
      for( qk in which(NonZeroCols)){
        temp.Trans[,qk]=temp.Trans[,qk] / sum(temp.Trans[,qk])
      }

    } else{
      qi = qxk[j - 1]
      qi.x = xk[, j - 1]

      for (qj in 1:N) {
        qj.x = ItoB(qj,n,K)
        # if (sum(qj.x) != sum(x[, j]))
        #   next
        if( !CheckMarginMatch(qj.x,x[, j],K) )
          next
        if (!CheckMatch(x[, j], x[, j - 1], qj.x,  qi.x,K) |
            !CheckMatch(x[, j], x[, j - 1], qj.x,  x[, j - 1],K))
          next

        if (j < J) {
          for (qk in 1:N) {
            qk.x = ItoB(qk,n,K)
            # if (sum(qk.x) != sum(x[, j + 1]))
            #   next
            if( !CheckMarginMatch(qk.x,x[, j + 1],K) )
              next
            if (!CheckMatch(x[, j], x[, j + 1], qj.x,  qk.x,K))
              next
            Debug.Inner(n)
            temp.Trans[qj, qk] = 1
            NonZeroCols[qk]=T
          }
        } else{
          # for j==J, no qk
          Debug.Inner(n)
          temp.Trans[qj, 1] = 1
          NonZeroCols[1]=T
        }

      }

      # inefficient
      # temp.Trans = apply(temp.Trans  , 2, function(v) {  # inefficient!!
      #   if (sum(v) > 0) {
      #     v =    lastTrans[qi ,] * v
      #     v / sum(v)
      #   }
      #   else
      #     0 * v
      # })
      for( qk in which(NonZeroCols)){
        v=temp.Trans[,qk] * lastTrans[ qi ,]
        temp.Trans[,qk]=v/sum(v)
      }

    }

    # draw a sample
    if (j == J) {
      qxk[j] = sample(N, 1, prob = temp.Trans[ , 1])
    } else{
      qxk[j] = sample(N, 1, prob = temp.Trans[ , qx[j + 1]])
    }
    lastTrans=temp.Trans
    xk[, j] = ItoB(qxk[j],n,K)

  }

  if(bMinus)xk[xk==0]=-1
  return(xk)
}


library(knockoff)

cknofkoff.MC.SCIP.split = function(x , n0=5, K=2 ,sparse=F) {
  n = nrow(x)
  J = ncol(x)
  xk=NULL

  AA=split(as.data.frame(x), rep(1:ceiling(n/n0), each=n0, length.out=n) ) # fast

  if(!sparse){
    xk=do.call(rbind,
               lapply(
                 AA,
                 cknofkoff.MC.SCIP, J=J,K=K)
    )
  }else{
    xk=do.call(rbind,
               lapply(
                 AA,
                 cknofkoff.MC.SCIP.sparse, J=J,K=K)
    )
  }

  xk
}
