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

fit = svd(ppmi)

saveRDS(fit,file=paste0('output/SVD.',source.name,'.rds'))
