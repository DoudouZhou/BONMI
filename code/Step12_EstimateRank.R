#Estimate of rank and weight
#How to choose rank: by the spearman similarity of a validation set
args <- commandArgs(TRUE)

s0 <- as.numeric(args[[1]])

if(s0==1){
  load('data/ppmi.stan.concepts_perBin_30d.RData')
  stan.dict = read.csv('data/singlets_concepts_perBin_30d.txt',
                       sep=',',header = T, stringsAsFactors = FALSE)
  cui.stan = stan.dict$cui[stan.dict$freq>1000]
  ppmi.stan = as.matrix(ppmi.stan)
  
  ppmi = ppmi.stan[stan.dict$freq>1000,stan.dict$freq>1000]
  
  source.name = 'stanford'
}
if(s0==2){
  load("data/biobank_sPMI_common.Rdata")
  ppmi = ppmi.bi
  ppmi = as.matrix(ppmi)
  source.name = 'biobank'
}
if(s0==3){
  load('data/ppmi.mimic.4.23.RData')
  ppmi = ppmi.mi
  ppmi = as.matrix(ppmi)
  source.name = 'mimic'
}
if(s0==4){
  Chpmi <- readRDS("data/PPMI_Other8_window_10threshold_0weighted_0.rds")
  ppmi = Chpmi
  ppmi = as.matrix(ppmi)
  source.name = 'chinese'
}

fit = readRDS(file=paste0('output/SVD.',source.name,'.rds'))

source('EvalFunction.R')

Res = NULL
for(r in seq(100,4000,by=100)){
  U = Em(fit,r)
  U = U/apply(U, 1, norm, '2')
  Cos = U%*%t(U)
  rownames(Cos) = colnames(Cos) = rownames(ppmi)
  sd = Sigma2(fit,r)
  
  res = c(r,sd)
  if(s0==4){
    res1 = CoSim(Cos)
    res2 = CoSim2(Cos)
    res = c(res,res1,res2)
  }else{
    res1 = CoSim_CUI(Cos)
    res2 = CoSim_CUI2(Cos)
    res3 = CoSim_CUI3(Cos)
    res = c(res,res1,res2,res3)
  }
  Res = rbind(Res,res)
}

saveRDS(Res,file=paste0('output/Step12_rank.',source.name,'.rds'))