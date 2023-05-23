#or we can input the rank

#Decide the weights
load('output/Step13_allvalidation.Rdata')
S1 = res_stanford$Sim1 + res_biobank$Sim1 + res_mimic$Sim1 + res_chinese$Sim1
print(which.max(S1))
r = res_stanford$rank[which.max(S1)]
print(r)
idr = which(res_stanford$rank==r)
weights = c(res_stanford$sd[idr] , res_biobank$sd[idr] , res_mimic$sd[idr] , res_chinese$sd[idr])
weights = 1/weights^2
print(round(weights,3))

#Load the PMI data
load('data/ppmi.stan.concepts_perBin_30d.RData')
stan.dict = read.csv('data/singlets_concepts_perBin_30d.txt',
                       sep=',',header = T, stringsAsFactors = FALSE)
cui.stan = stan.dict$cui[stan.dict$freq>1000]
ppmi.stan = as.matrix(ppmi.stan)
ppmi.stan = ppmi.stan[stan.dict$freq>1000,stan.dict$freq>1000]
  
load("data/biobank_sPMI_common.Rdata")
ppmi.bi = as.matrix(ppmi.bi)
cui.bio = rownames(ppmi.bi)

load('data/ppmi.mimic.4.23.RData')
ppmi.mi = as.matrix(ppmi.mi)
cui.mimic = rownames(ppmi.mi)

Chpmi <- readRDS("data/PPMI_Other8_window_10threshold_0weighted_0.rds")
ppmi.ch = as.matrix(Chpmi)
cui.ch0 = rownames(ppmi.ch)
#########
load('data/TranslationLabel.RData')
cui.ch = cui.ch0; cui.ch[match(Train$Chinese,cui.ch0)] = Train$CUI
rownames(ppmi.ch) = cui.ch

Cuiall = unique(c(cui.stan, cui.bio, cui.mimic, cui.ch))
N = length(Cuiall)

Cuiall.col = Cuiall
Cuiall.col[match(Train$CUI,Cuiall)] = Train$Chinese

save(cui.stan, cui.bio, cui.mimic,cui.ch0,cui.ch, Cuiall,Cuiall.col,
     file='output/Step21_RownamesAndColanmes.Rdata')

###Combine the data
W = list()
W[[1]] = ppmi.stan; W[[2]] = ppmi.bi; W[[3]] = ppmi.mi; W[[4]] = ppmi.ch
saveRDS(W,file='output/Step21_W.Rds')

source('Gdata-MS.R')
source('functions-MS.R')

#Ww includes the integrated big W (with missing entries completed by zero (Wc); or NA (Wo) )
#and the weight matrix Weis and Omega (C)
Ww = CreateW(W,weights,Cuiall)

datnew = New_W(W,Ww$Wc)

save(Ww,datnew,file='output/Step21_Aggregation.Rdata')

print(table(Ww$Weis))
print(mean(C))


####Formal Version
load('output/Step13_allvalidation.Rdata')
load(file='output/Step21_RownamesAndColanmes.Rdata')

W = readRDS(file='output/Step21_W.Rds')

source('Gdata-MS.R')
source('functions-MS.R')

#c(2200,2900,3500)
for(r in 300){
  print(r)
  idr = which(res_stanford$rank==r)
  weights = c(res_stanford$sd[idr] , res_biobank$sd[idr] , res_mimic$sd[idr] , res_chinese$sd[idr])
  weights = 1/weights^2
  print(round(weights,3))
  
  Ww = CreateW(W,weights,Cuiall)
  datnew = New_W(W,Ww$Wc)
  
  save(Ww,datnew,file=paste0('output/Step21_Aggregation_',r,'.Rdata'))
}
