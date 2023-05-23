##May
# source('Step21_Aggregation.R')

args <- commandArgs(TRUE)

r <- as.numeric(args[[1]])
meth <- as.numeric(args[[2]])
lambda <- as.numeric(args[[3]])
#r = 100

source('Gdata-MS.R')
source('functions-MS.R')
#load('output/Step21_Aggregation.Rdata')
load(paste0('output/Step21_Aggregation_',r,'.Rdata'))
print(str(datnew))
print(meth)

if(meth==1){
  mod.Orth <- Orth(datnew$W,Ww$Wc,Ww$C,r)
  saveRDS(mod.Orth,file = paste0('output/Step22_BELT_r',r,'.Rds'))
}

if(meth==2){
  mod.zero <- Zero(Ww$Wc,r)
  saveRDS(mod.zero,file = paste0('output/Step22_Zero_r',r,'.Rds'))
}

if(meth==3){
  mod.Anru <- SMC(datnew$W,Ww$Wc,Ww$C,r)
  saveRDS(mod.Anru,file = paste0('output/Step22_SMC_r',r,'.Rds'))
}

if(meth==4){
  p = sum(Ww$C>0)/nrow(Ww$C)^2
  # fit.W = readRDS(file = paste0('output/Step22_Zero_r',r,'.Rds'))
  # fit.W = fit.W$fit
  set.seed(1)
  fit.W = rsvd(Ww$Wc,r)
  tau = fit.W$d[1]/fit.W$d[r]
  sigma_max = fit.W$d[1]/p
  
  #stepsize = 2/(25*tau*sigma_max)
  stepsize = 0.1/(25*tau*sigma_max)
  X = fit.W$u[,1:r]%*%diag(sqrt(fit.W$d[1:r]))
  
  mod.Vgd = Vgd(Ww$Wc,Ww$C,X0=X,stepsize,Tmax=100)
  saveRDS(mod.Vgd,file = paste0('output/Step22_Vgd_r',r,'.Rds'))
}

if(meth==5){
  # fit.W = readRDS(file = paste0('output/Step22_Zero_r',r,'.Rds'))
  # fit.W = fit.W$fit
  # mod.ALS = softImpute(Ww$Wo, rank.max = r, lambda = lambda, type = "als", thresh = 1e-05,
  #                      maxit = 50, trace.it = FALSE, warm.start = fit.W, final.svd = TRUE)
  
  mod.ALS = softImpute(Ww$Wo, rank.max = r, lambda = lambda, type = "als", thresh = 1e-05,
                       maxit = 50, trace.it = FALSE, final.svd = TRUE)
  print(length(mod.ALS$d))
  saveRDS(mod.ALS,file = paste0('output/Step22_ALS_r',r,'_lambda',lambda,'.Rds'))
}


if(meth==6){
  p = sum(Ww$C>0)/nrow(Ww$C)^2
  # fit.W = readRDS(file = paste0('output/Step22_Zero_r',r,'.Rds'))
  # fit.W = fit.W$fit
  set.seed(1)
  fit.W = readRDS(paste0('output/Step22_BELT_r',r,'.Rds'))
  tau = fit.W$fit$d[1]/fit.W$fit$d[r]
  sigma_max = fit.W$fit$d[1]
  
  #stepsize = 2/(25*tau*sigma_max)
  stepsize = 0.05/(25*tau*sigma_max)
  X = fit.W$X
  
  Vgd0<-function(W,C,X0,stepsize,Tmax){
    t=0
    Omega = C
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
  
  
  Ww$C[Ww$C==0] = 0.5
  mod.Vgd = Vgd0(Ww$Wc,Ww$C,X0=X,stepsize,Tmax=100)
  saveRDS(mod.Vgd,file = paste0('output/Step22_Vgd0_r',r,'.Rds'))
}

if(meth==7){
  print(meth)
  fit.W = readRDS(file = paste0('output/Step22_BELT_r',r,'.Rds'))
  fit.W = fit.W$fit
  print('Begin')
  mod.ALS = softImpute(Ww$Wo, rank.max = r, lambda = lambda, type = "als", thresh = 1e-05,
                       maxit = 50, trace.it = TRUE, warm.start = fit.W, final.svd = TRUE)
  print(length(mod.ALS$d))
  saveRDS(mod.ALS,file = paste0('output/Step22_BELTandALS_r',r,'_lambda',lambda,'.Rds'))
  
}