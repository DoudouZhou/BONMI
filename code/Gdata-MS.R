##Create a big matrix arranged by namestar
#W is the list of the m observed matrices
CreateW<-function(W,weights,namestar){
  m = length(W)
  
  Wc = matrix(0,nrow=length(namestar),ncol=length(namestar))
  C = matrix(0,nrow=length(namestar),ncol=length(namestar))
  Weis = matrix(0,nrow=length(namestar),ncol=length(namestar))
  
  rownames(Wc) = colnames(Wc) = namestar
  rownames(C) = colnames(C) = namestar
  rownames(Weis) = colnames(Weis) = namestar
  for(s in 1:m){
    id = match(rownames(W[[s]]),namestar)
    Wc[id,id] = Wc[id,id] + weights[s]*W[[s]]
    C[id,id] = C[id,id]+1
    Weis[id,id] = Weis[id,id]+weights[s] 
  }
  Wc[C>0] = Wc[C>0]/Weis[C>0]
  Wo = Wc; Wo[C==0] = NA
  return(list('Wc'=Wc,'Wo'=Wo,'C'=C,'Weis'=Weis))
}

New_W<-function(W,Wc){
  m = length(W)
  W.new = list()
  for(s in 1:m){
    Is = match(rownames(W[[s]]),rownames(Wc))
    W.new[[s]] = Wc[Is,Is]
  }
  return(list(W=W.new))
}

Obtain_weight<-function(W,r){
  SVDs = lapply(1:m, function(s){svd(W[[s]],nu=r,nv=r)})
  Sigmas = sapply(1:m, function(s){
   norm(W[[s]]-SVDs[[s]]$u[,1:r]%*%diag(SVDs[[s]]$d[1:r])%*%t(SVDs[[s]]$u[,1:r]),'F')/nrow(W[[s]])})
  1/Sigmas^2
}
