library(softImpute)
library(rsvd)
#library(rsparse)

#Learn the orthogonal transformation matrix O
Procrustes <- function(X1,X2){
  #Return Omega = arg min||X1 - X2 Omgea||_F
  H = t(X2)%*%X1
  mod = svd(H)
  return(mod$u%*%t(mod$v))
}

embedding<-function(svd,r){
  svd$u[,1:r]%*%diag(sqrt(svd$d[1:r]))
}

Orth<-function(W,Wc,C,r){
  #C[i,j]=0 means the (i,j) entry is missed
  m = length(W)
  SVDs = lapply(1:m, function(s){svd(W[[s]],nu=r,nv=r)})
  Xs = lapply(1:m, function(s){ 
    U = embedding(SVDs[[s]],r)
    rownames(U) = rownames(W[[s]])
    U
  })
  Wm = matrix(0,nrow=nrow(Wc),ncol=ncol(Wc))
  M = matrix(0,nrow=nrow(Wc),ncol=ncol(Wc))
  rownames(Wm) = colnames(Wm) = rownames(Wc) 
  rownames(M) = colnames(M) = rownames(Wc) 
  
  for(s in 1:(m-1)){
    for(k in (s+1):m){
      name12 = intersect(rownames(W[[s]]),rownames(W[[k]]))
      ids = match(name12,rownames(W[[s]]))
      idk = match(name12,rownames(W[[k]]))
      Osk = Procrustes(Xs[[s]][ids,],Xs[[k]][idk,])
      Wsk = Xs[[s]][-ids,]%*%t(Osk)%*%t(Xs[[k]][-idk,])
      id1 = match(rownames(W[[s]])[-ids],rownames(Wm))
      id2 = match(rownames(W[[k]])[-idk],rownames(Wm))
      Wm[id1,id2] = Wm[id1,id2]+Wsk; Wm[id2,id1] =  Wm[id2,id1]+t(Wsk)
      M[id1,id2] = M[id1,id2] + 1; M[id2,id1] = M[id2,id1] + 1
    }
  }
  Wm[M>0] = Wm[M>0]/M[M>0]
  Wm[C>0] = Wc[C>0]
  M[C>0] = 0
  
  #fit.W = svd(Wm,nu=r,nv=r)
  set.seed(1)
  fit.W = rsvd(Wm,r+2000)
  X = embedding(fit.W,r)
  What = X%*%t(X)
  rownames(X) = rownames(What) = colnames(What) = rownames(Wm)
  return(list('Wm'=Wm,'M'=M,
              'X'=X,'What'=What,'fit'=fit.W,'Xs'=Xs))
}

smc<-function(A11,A12,A21,r){
  p1 = nrow(A11)+nrow(A21)
  m1 = nrow(A11)
  
  Adot1 = rbind(A11,A21)
  A1dot = cbind(A11,A12)
  
  #
  set.seed(1)
  fit1 = rsvd(Adot1,r)
  fit2 = rsvd(A1dot,r)
  
  Z11 = t(fit2$u)%*%A11%*%fit1$v
  Z12 = t(fit2$u)%*%A12
  Z21 = A21%*%fit1$v

  A22 = Z21[,1:r]%*%solve(Z11[1:r,1:r])%*%Z12[1:r,]
  rownames(A22) = rownames(A21); colnames(A22) = colnames(A12)
  A22
}

SMC<-function(W,Wc,C,r){
  #C[i,j]=0 means the (i,j) entry is missed
  m = length(W)
  
  Wm = matrix(0,nrow=nrow(Wc),ncol=ncol(Wc))
  M = matrix(0,nrow=nrow(Wc),ncol=ncol(Wc))
  rownames(Wm) = colnames(Wm) = rownames(Wc) 
  rownames(M) = colnames(M) = rownames(Wc) 
  
  for(s in 1:(m-1)){
    for(k in (s+1):m){
      name12 = intersect(rownames(W[[s]]),rownames(W[[k]]))
      ids = match(name12,rownames(W[[s]]))
      idk = match(name12,rownames(W[[k]]))
      W1 = W[[s]]; W2 = W[[k]]
      
      A11 = W2[,idk]; A11[idk,] = (A11[idk,]+W1[ids,ids])/2
      A12 = W2[,-idk]; A21 = W1[-ids,ids]
      A22.r = smc(A11,A12,A21,r)
      
      A11 = W1[,ids]; A11[ids,] = (A11[ids,]+W2[idk,idk])/2
      A12 = W1[,-ids]; A21 = W2[-idk,idk]
      A22.l = smc(A11,A12,A21,r)
      
      A22 = (A22.r + t(A22.l))/2
      
      id1 = match(rownames(W[[s]])[-ids],rownames(Wm))
      id2 = match(rownames(W[[k]])[-idk],rownames(Wm))
      Wm[id1,id2] = Wm[id1,id2]+A22; Wm[id2,id1] =  Wm[id2,id1]+t(A22)
      M[id1,id2] = M[id1,id2] + 1; M[id2,id1] = M[id2,id1] + 1
    }
  }
  Wm[M>0] = Wm[M>0]/M[M>0]
  Wm[C>0] = Wc[C>0]
  M[C>0] = 0
  fit.W = rsvd(Wm,r+100)
  X = embedding(fit.W,r)
  What = X%*%t(X)
  rownames(X) = rownames(What) = colnames(What) = rownames(Wm)
  return(list('Wm'=Wm,'M'=M,
              'X'=X,'What'=What,'fit'=fit.W))
}

Zero<-function(Wc,r){
  set.seed(1)
  #fit.W = svd(Wc,nu=r,nv=r)
  fit.W = rsvd(Wc,r+2000)
  X = embedding(fit.W,r)
  What = X%*%t(X)
  rownames(X) = rownames(What) = rownames(Wc)
  return(list('X'=X,'What'=What,'fit'=fit.W))
}

Vgd<-function(W,C,X0,stepsize,Tmax){
  p = sum(C>0)/nrow(C)^2
  
  #Spectral initialization
  #X0 = X0/sqrt(p)
  #Gradient updates
  t=0
  Omega = 1*(C>0)
  repeat({
    t = t+1
    if(t>Tmax){
      break
    }
    grad = (X0%*%t(X0)*Omega - W)%*%X0/p
    X0 = X0 - stepsize*grad
    relerr = norm(grad,'F')/norm(X0,'F')
    print(c(t,round(relerr,3)))
    if(relerr<1e-2){
      break
    }
  })
  What = X0%*%t(X0)
  return(list('What'=What,'X'=X0,'stepsize'=stepsize))
}

#Obtain Xs for the original machine translation method
Original<-function(Xs){
  m = length(Xs)
  for(s in 2:m){
    train = intersect(rownames(Xs[[1]]),rownames(Xs[[s]]))
    id1 = match(train,rownames(Xs[[1]]))
    ids = match(train,rownames(Xs[[s]]))
    O = Procrustes(Xs[[1]][id1,],Xs[[s]][ids,])
    Xs[[s]] = Xs[[s]]%*%O
  }
  Xs
}

#W = list(W1,...,Wm)
BuildXs<-function(W,X){
  m = length(W)
  Xs = lapply(1:m,function(s){
    X[match(rownames(W[[s]]),rownames(X)),]
  })
  Xs
}
#Input: m sets of embeddings: Xs^{test} \in R^{(ns + ntest) x r}
library(stringr)
Translation<-function(Xs){
  m = length(Xs)
  Xs[[1]] = Xs[[1]]/apply(Xs[[1]],1, norm, '2')
  id1 = which(!is.na(str_match(rownames(Xs[[1]]),'Test')[,1]) )
  Acc = NULL
  for(s in 2:m){
    Xs[[s]] = Xs[[s]]/apply(Xs[[s]],1, norm, '2')
    ids = which(!is.na(str_match(rownames(Xs[[1]]),'Test')[,1]) )
    costest = Xs[[s]][ids,]%*%t(Xs[[1]])
    map = t(apply(costest, 1, order,decreasing = TRUE))
    acc = sapply(c(1,3,5), function(k){
      sum(sapply(1:nrow(costest), function(i){ids[i]%in%map[i,1:k]}))
    })/nrow(costest)
    Acc = rbind(Acc,acc)
  }
  colMeans(Acc)
}
