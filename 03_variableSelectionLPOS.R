# Name: 03_variableSelectionLPOS.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 27/02/2018
# Desc: variable selection using the mass spec data

source('header.R')

if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/experimental/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

dfData = read.csv(file.choose(), header = T, na.strings = c('na', 'NA', 'NaN'))
dfData = t(dfData)
dfData = data.frame(dfData)
cn = t(dfData[1,])
colnames(dfData) = cn[,1]
dfData = dfData[-1,]
rm(cn)

## remove the unwanted groups 
i = which(dfData$Class_ID %in% c('ACLF', 'AD', 'SC'))
length(i)
dfData = dfData[i,]
dim(dfData)
dfData = droplevels.data.frame(dfData)

gc(reset = T)

write.csv(dfData, 'temp/transposed_cleaned_lpos.csv', row.names = F)
## clears data types etc, do a fresh reload
dfData = read.csv('temp/transposed_cleaned_lpos.csv', header = T, na.strings = c('na', 'NA', 'NaN'))
gc(reset = T)
## main grouping factor
fGroups = factor(dfData$X90_day)
levels(fGroups)

colnames(dfData)[1:100]
## first 94 columns are not mass spec data
## drop those first to work only on massspec data
dfData = dfData[,-c(1:94)]
gc(reset = T)
dim(dfData)
dfData = na.omit(dfData)
dim(dfData)

## log transform the data
mDat = log(dfData+1)
mDat = t(mDat)

#### calculate a scaling factor for normalization
## this outlier detection is based on previous analysis
m = colMeans(mDat)
plot(m)
## 3 samples are different or outliers drop those
k = hclust(dist(m))
plot(k)
c = cutree(k, k = 2)
table(c)
#iOutliers = which(c == 2)

## drop the outliers
#mDat = mDat[,-iOutliers]
sf = rowMeans(mDat)
mDat.res = sweep(mDat, 1, sf, '-')
## use median as size factor
sf = apply(mDat.res, 2, function(x) quantile(x, prob=0.5))
mDat.norm = sweep(mDat, 2, sf, '-')

## this comparison of normalised and raw data was not done for this data set
## compare the normalised and raw data
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

oDiag.1 = CDiagnosticPlots(mDat.norm, 'Normalised')
oDiag.2 = CDiagnosticPlots(mDat, 'Original')
# the batch variable we wish to colour by, 
# this can be any grouping/clustering in the data capture process
# e.g. in this case it is different lanes/machines
fBatch = fGroups

## compare the 2 methods using various plots
par(mfrow=c(1,2))
boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)
boxplot.median.summary(oDiag.2, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)

plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.mean.summary(oDiag.2, fBatch, axis.label.cex = 0.7)

plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.sigma.summary(oDiag.2, fBatch, axis.label.cex = 0.7)

#plot.missing.summary(oDiag.1, fBatch, axis.label.cex = 0.7, cex.main=1)

# plot.PCA(oDiag.1, fBatch, cex.main=1)
# 
# plot.dendogram(oDiag.1, fBatch, labels_cex = 0.8, cex.main=0.7)

## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag.1)
l$PCA.jitter = F
l$HC.jitter = F

oDiag.1 = CDiagnosticPlotsSetParameters(oDiag.1, l)
oDiag.2 = CDiagnosticPlotsSetParameters(oDiag.2, l)
plot.PCA(oDiag.1, fBatch)
plot.PCA(oDiag.2, fBatch)
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.7)
plot.dendogram(oDiag.2, fBatch, labels_cex = 0.7)

## comparisons of normalised and raw data done

## perform analysis with normalised
mDat = mDat.norm
dim(mDat)
rm(mDat.res); rm(mDat.norm)

## reload the full data
dfData = read.csv('temp/transposed_cleaned_lpos.csv', header = T, na.strings = c('na', 'NA', 'NaN'))
dfData = dfData[,-c(1,2)]
dim(dfData)
## add the new normalised/raw data
colnames(dfData)[1:100]
dfData = dfData[,1:92]
dim(mDat)
dim(dfData)
identical(rownames(dfData), colnames(mDat))
dfData = cbind(dfData, t(mDat))
dim(dfData)

gc(reset = T)

fGroups = factor(dfData$X90_day)
## control variable
fControl = factor(dfData$Class_ID)
dfData = dfData[,-(which(colnames(dfData) %in% c('Class_ID', 'X90_day', 'X30_day', 'X1_Year', "Sample.File.Name")))]

# save backup
dfData.bk = dfData
fGroups.bk = fGroups
fControl.bk = fControl
## perform variable selection only on the mass spec data
colnames(dfData)[1:100]
dfData = dfData[,-(1:88)]
dim(dfData)

## create a test and training set
set.seed(1234);
test = sample(1:nrow(dfData), size = nrow(dfData)*0.2, replace = F)
table(fGroups[test])
table(fGroups[-test])

table(fControl[test])
table(fControl[-test])

### DE model using limma
library(limma)

design = model.matrix(~ fGroups[-test] + fControl[-test])
head(design)

mData = as.matrix(dfData[-test,])

fit = lmFit(t(mData), design)
fit = eBayes(fit)

dfLimmma.2 = topTable(fit, coef = 2, adjust='BH', number=Inf)
hist(dfLimmma.2$logFC)
hist(dfLimmma.2$adj.P.Val)
dfLimmma.2 = dfLimmma.2[order(dfLimmma.2$P.Value, decreasing = F),]
head(dfLimmma.2)

########## repeat this analysis in stan using t distribution hierarchical model
## format the data frame for input
dfData.org = dfData
dfData = stack(dfData[-test, ])
dfData$fBatch = fGroups[-test]
dfData$fAdjust = fControl[-test]
dfData$Coef = factor(dfData$fBatch:dfData$ind)
dfData$Coef.adj = factor(dfData$fAdjust:dfData$ind)
dfData = droplevels.data.frame(dfData)
dfData = dfData[order(dfData$Coef, dfData$Coef.adj), ]

# #### fit mixed effect model
library(lme4)
fit.lme1 = lmer(values ~ 1 + (1 | Coef) + (1 | Coef.adj), data=dfData, REML=F)
summary(fit.lme1)

# plot(fitted(fit.lme1), resid(fit.lme1), pch=20, cex=0.7)
# lines(lowess(fitted(fit.lme1), resid(fit.lme1)), col=2)

## fit model with stan
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

################# t model
stanDso = rstan::stan_model(file='tResponse2RandomEffectNoFixed.stan')

## calculate hyperparameters for variance of coefficients
l = gammaShRaFromModeSD(sd(dfData$values), 2*sd(dfData$values))

## set initial values
# initf = function(chain_id = 1) {
#   gm = tapply(dfData$values, dfData$Coef, mean) - mean(dfData$values)
#   list(betas = mean(dfData$values), sigmaRan1 = sd(gm), sigmaPop=sd(dfData$values), nu=4, rGroupsJitter1=gm)
# }
## set initial values
ran = ranef(fit.lme1)
r1 = ran$Coef
r2 = ran$Coef.adj
initf = function(chain_id = 1) {
  list(sigmaRan1 = 2, sigmaRan2=2, sigmaPop=1, rGroupsJitter1=r1, rGroupsJitter2=r2, nu=4)
}

### try a t model without mixture
lStanData = list(Ntotal=nrow(dfData), Nclusters1=nlevels(dfData$Coef),
                 Nclusters2=nlevels(dfData$Coef.adj),
                 NgroupMap1=as.numeric(dfData$Coef),
                 NgroupMap2=as.numeric(dfData$Coef.adj),
                 Ncol=1, 
                 y=dfData$values, 
                 gammaShape=l$shape, gammaRate=l$rate,
                 intercept = mean(dfData$values), intercept_sd= sd(dfData$values)*3)

fit.stan = sampling(stanDso, data=lStanData, iter=500, chains=2,
                    pars=c('betas', 'sigmaRan1', 'sigmaRan2',
                           'nu', 'sigmaPop', #'mu',
                           'rGroupsJitter1', 'rGroupsJitter2'),
                    cores=2, init=initf)#, control=list(adapt_delta=0.99, max_treedepth = 12))
save(fit.stan, file='temp/fit.stan.tdis_lpos_27_Feb_500.rds')
print(fit.stan, c('betas', 'sigmaRan1', 'sigmaRan2', 'sigmaPop', 'nu'), digits=3)

### get the coefficient for main treatment 
mCoef = extract(fit.stan)$rGroupsJitter1
dim(mCoef)
# ## get the intercept at population level
iIntercept = as.numeric(extract(fit.stan)$betas)
## add the intercept to each random effect variable, to get the full coefficient
mCoef = sweep(mCoef, 1, iIntercept, '+')

## function to calculate statistics for differences between coefficients
getDifference = function(ivData, ivBaseline){
  stopifnot(length(ivData) == length(ivBaseline))
  # get the difference vector
  d = ivData - ivBaseline
  # get the z value
  z = mean(d)/sd(d)
  # get 2 sided p-value
  p = pnorm(-abs(mean(d)/sd(d)))*2
  return(list(z=z, p=p))
}

## split the data into the comparisons required
d = data.frame(cols=1:ncol(mCoef), mods=levels(dfData$Coef))
## split this factor into sub factors
f = strsplit(as.character(d$mods), ':')
d = cbind(d, do.call(rbind, f))
head(d)
colnames(d) = c(colnames(d)[1:2], c('fBatch', 'ind'))
d$split = factor(d$ind)

## get a p-value for each comparison
l = tapply(d$cols, d$split, FUN = function(x, base='0', deflection='1') {
  c = x
  names(c) = as.character(d$fBatch[c])
  dif = getDifference(ivData = mCoef[,c[deflection]], ivBaseline = mCoef[,c[base]])
  r = data.frame(ind= as.character(d$ind[c[base]]), coef.base=mean(mCoef[,c[base]]), 
                 coef.deflection=mean(mCoef[,c[deflection]]), zscore=dif$z, pvalue=dif$p)
  r$difference = r$coef.deflection - r$coef.base
  #return(format(r, digi=3))
  return(r)
})

dfResults = do.call(rbind, l)
dfResults$adj.P.Val = p.adjust(dfResults$pvalue, method='BH')

### compare the results from the 2 models
dfResults$logFC = dfResults$difference
dfResults$P.Value = dfResults$pvalue
dfLimmma.2$SYMBOL = as.character(rownames(dfLimmma.2))
dfResults$SYMBOL = as.character(rownames(dfResults))

## produce the plots 
f_plotVolcano(dfLimmma.2, 'limma 1 vs 0', fc.lim = c(-2, 2), p.adj.cut = 1)
f_plotVolcano(dfResults, 'Stan 1 vs 0')

m = tapply(dfData$values, dfData$ind, mean)
i = match(rownames(dfResults), names(m))
m = m[i]
identical(names(m), rownames(dfResults))
plotMeanFC(m, dfResults, 0.1, 'Stan 1 vs 0')

m = tapply(dfData$values, dfData$ind, mean)
i = match(rownames(dfLimmma.2), names(m))
m = m[i]
m = m[!is.na(m)]
i = match(names(m), rownames(dfLimmma.2))
dfLimmma.2 = dfLimmma.2[i,]
identical(names(m), rownames(dfLimmma.2))
plotMeanFC(m, dfLimmma.2, 0.1, 'limma 1 vs 0')

i = match(rownames(dfResults), rownames(dfLimmma.2))
dfLimmma.2 = dfLimmma.2[i,]
i = match(rownames(dfLimmma.2), rownames(dfResults))
dfResults = dfResults[i,]
identical(rownames(dfResults), rownames(dfLimmma.2))

plot(dfResults$pvalue, dfLimmma.2$P.Value, pch=20, cex=0.6, col='grey', main='P Values 1 vs 0', xlab='Stan', ylab='Limma')
abline(lm(dfLimmma.2$P.Value ~ dfResults$pvalue), col=2, lwd=2)
plot(dfResults$logFC, dfLimmma.2$logFC, pch=20, cex=0.8, col='grey', main='Log FC 1 vs 0', xlab='Stan', ylab='Limma')
abline(lm(dfLimmma.2$logFC ~ dfResults$logFC), col=2, lwd=1)
df = cbind(stan=dfResults$pvalue, limma=dfLimmma.2$P.Value)

write.csv(dfResults, file='temp/stan_t_lpos.csv', row.names = F)
write.csv(dfLimmma.2, file='temp/limma_lpos.csv', row.names = F)

######### stan section ends
dfData = dfData.org
dfResults = dfResults[order(dfResults$pvalue), ]
table(dfResults$adj.P.Val < 0.1)
## select the top variables 1 to 2000 ranked by pvalue
cvTopVariables = rownames(dfResults)[1:2000]#[dfResults$adj.P.Val < 0.1]
length(cvTopVariables)
## perform nested random forest
## adjust boot.num as desired
# oVar.r = CVariableSelection.RandomForest(dfData[-test, cvTopVariables], fGroups[-test], boot.num = 100, big.warn = F)
# save(oVar.r, file='temp/oVar_lpos.rds')
# # plot the top 20 variables based on importance scort with 95% confidence interval for standard error
# par(mfrow=c(1,1))
# plot.var.selection(oVar.r)
# # get the variables
# dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
# # select the top 30 variables
# cvTopGenes = rownames(dfRF)[1:30]

# use the top 30 features to find top combinations of genes
# dfData = dfData[,colnames(dfData) %in% cvTopGenes]

# ## look at colinear variables
# m = NULL;
# 
# for (i in 1:ncol(dfData)){
#   m = cbind(m, dfData[-test ,i])
# }
# colnames(m) = colnames(dfData)
# mCor = cor(m, use="na.or.complete")
# library(caret)
# ### find the columns that are correlated and should be removed
# n = findCorrelation((mCor), cutoff = 0.7, names=T)
# data.frame(n)
# sapply(n, function(x) {
#   (abs(mCor[,x]) >= 0.7)
# })
# s = sapply(n, function(x) {
#   (abs(mCor[,x]) >= 0.7)
# })
# colSums(s)
# cvKeep = names(colSums(s)[colSums(s) <= 4])
# #cvKeep = c('SOFA', 'MELD', 'Albumin', 'BMI', 'Neutrophil')
# n = n[!(n%in% cvKeep)]
# i = which(colnames(dfData) %in% n)
# cn = colnames(dfData)[-i]
# 
# dfData.bk2 = dfData
# dfData = dfData[,cn]
# dim(dfData)

# oVar.sub = CVariableSelection.ReduceModel(dfData[-test, ], fGroups[-test], boot.num = 100)
# # plot the number of variables vs average error rate
# plot.var.selection(oVar.sub)
# 
# # print variable combinations
# for (i in 1:6){
#   cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
#   cat('Variable Count', i, paste(cvTopGenes.sub), '\n')
#   #print(cvTopGenes.sub)
# }
# 
# ## 10 fold nested cross validation with various variable combinations
# par(mfrow=c(2,2))
# # try models of various sizes with CV
# for (i in 1:6){
#   cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
#   dfData.train = data.frame(dfData[-test ,cvTopGenes.sub])
#   colnames(dfData.train) = cvTopGenes.sub
#   
#   dfData.test = data.frame(dfData[test ,cvTopGenes.sub])
#   colnames(dfData.test) = cvTopGenes.sub
#   
#   oCV = CCrossValidation.LDA(test.dat = dfData.test, train.dat = dfData.train, test.groups = fGroups[test],
#                              train.groups = fGroups[-test], level.predict = '0', boot.num = 500)
#   
#   plot.cv.performance(oCV)
#   # print variable names and 95% confidence interval for AUC
#   temp = oCV@oAuc.cv
#   x = as.numeric(temp@y.values)
#   print(paste('Variable Count', i))
#   print(cvTopGenes.sub)
#   print(signif(quantile(x, probs = c(0.025, 0.975)), 2))
# }
##################################

############## use a binomial regression approach this time to rank variables
dfData = dfData.org
dfData = dfData[,cvTopVariables]
dfData = scale(dfData)
dim(dfData)
dfData = data.frame(dfData, fGroups)
dfData = dfData[-test, ]
#dfData$fGroups = fGroups

lData = list(resp=ifelse(dfData$fGroups == '0', 0, 1), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='binomialRegression.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

## give initial values
initf = function(chain_id = 1) {
  list(betas=rep(0, times=ncol(lStanData$X)), tau=0.5)
}


fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=2, pars=c('tau', 'betas2'), init=initf, 
                    control=list(adapt_delta=0.99, max_treedepth = 11))

print(fit.stan, c('betas2', 'tau'))
print(fit.stan, 'tau')
traceplot(fit.stan, 'tau')

## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$betas2
dim(mCoef)
# ## get the intercept at population level
iIntercept = mCoef[,1]
mCoef = mCoef[,-1]
# ## add the intercept to each coefficient, to get the full coefficient
# mCoef = sweep(mCoef, 1, iIntercept, '+')

## function to calculate statistics for a coefficient
getDifference = function(ivData){
  # get the difference vector
  d = ivData
  # get the z value
  z = mean(d)/sd(d)
  # get 2 sided p-value
  p = pnorm(-abs(mean(d)/sd(d)))*2
  return(p)
}

ivPval = apply(mCoef, 2, getDifference)
hist(ivPval)
plot(colMeans(mCoef), ivPval, pch=19)
m = colMeans(mCoef)
names(m) = colnames(lData$mModMatrix)[2:length(m)]
m = abs(m)
m = sort(m, decreasing = T)
#i = which(ivPval < 0.85)
# colnames(lData$mModMatrix)[i+1]
# colMeans(mCoef)[i]
#cvTopGenes.binomial = colnames(lData$mModMatrix)[i+1]
cvTopGenes.binomial = names(m)[1:30] #names(m[m > 0.25])

# ## where are these in the random forest table
# which(rownames(dfRF) %in% cvTopGenes.binomial)
# 
# cvTopGenes.first = NULL;
# for (i in 1:5){
#   cvTopGenes.first = append(cvTopGenes.first, CVariableSelection.ReduceModel.getMinModel(oVar.sub, i))
# }
# cvTopGenes.first = unique(cvTopGenes.first)
# length(cvTopGenes.first)
################### repeat the selection process again 
# select the top 30 variables
#cvTopGenes = rownames(dfRF)[1:30]
cvTopGenes = cvTopGenes.binomial; #unique(c(cvTopGenes.first, cvTopGenes.binomial))
#cvTopGenes = cvTopGenes.binomial
length(cvTopGenes)
#cvTopGenes = cvTopGenes[1:30]
# use the top 30 features to find top combinations of genes
dfData = dfData.org
dfData = dfData[,colnames(dfData) %in% cvTopGenes]
dim(dfData)
# ## look at colinear variables
# m = NULL;
# 
# for (i in 1:ncol(dfData)){
#   m = cbind(m, dfData[-test ,i])
# }
# colnames(m) = colnames(dfData)
# mCor = cor(m, use="na.or.complete")
# library(caret)
# ### find the columns that are correlated and should be removed
# n = findCorrelation((mCor), cutoff = 0.7, names=T)
# data.frame(n)
# sapply(n, function(x) {
#   (abs(mCor[,x]) >= 0.7)
# })
# s = sapply(n, function(x) {
#   (abs(mCor[,x]) >= 0.7)
# })
# colSums(s)
# cvKeep = names(colSums(s)[colSums(s) <= 4])
# n = n[!(n%in% cvKeep)]
# i = which(colnames(dfData) %in% n)
# cn = colnames(dfData)[-i]
# 
# dfData.bk2 = dfData
# dfData = dfData[,cn]
# dim(dfData)
oVar.sub.first = oVar.sub

oVar.sub = CVariableSelection.ReduceModel(dfData[-test, ], fGroups[-test], boot.num = 30)
# plot the number of variables vs average error rate
plot.var.selection(oVar.sub)

# print variable combinations
for (i in 1:7){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  cat('Variable Count', i, paste(cvTopGenes.sub), '\n')
  #print(cvTopGenes.sub)
}

## 10 fold nested cross validation with various variable combinations
par(mfrow=c(2,2))
# try models of various sizes with CV
for (i in 1:8){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  dfData.train = data.frame(dfData[-test ,cvTopGenes.sub])
  colnames(dfData.train) = cvTopGenes.sub
  
  dfData.test = data.frame(dfData[test ,cvTopGenes.sub])
  colnames(dfData.test) = cvTopGenes.sub
  
  oCV = CCrossValidation.LDA(test.dat = dfData.test, train.dat = dfData.train, test.groups = fGroups[test],
                             train.groups = fGroups[-test], level.predict = '0', boot.num = 500)
  
  plot.cv.performance(oCV)
  # print variable names and 95% confidence interval for AUC
  temp = oCV@oAuc.cv
  x = as.numeric(temp@y.values)
  print(paste('Variable Count', i))
  print(cvTopGenes.sub)
  print(signif(quantile(x, probs = c(0.025, 0.975)), 2))
}

### plot these genes of interest
dfData = dfData[,colnames(dfData) %in% CVariableSelection.ReduceModel.getMinModel(oVar.sub, 7)]
dim(dfData)

dfData = stack(dfData)
dfData$fBatch = fGroups
dfData$fAdjust = fControl
dfData$Coef = factor(dfData$fBatch:dfData$ind)
dfData$Coef.adj = factor(dfData$fAdjust:dfData$ind)
dfData = droplevels.data.frame(dfData)
str(dfData)

library(lattice)
densityplot(~ values | ind, groups=fBatch, data=dfData, auto.key = list(columns=2), pch=20, cex=0.5)
bwplot(values ~ fBatch | ind, data=dfData, type='b', panel=panel.violin, varwidth=F, xlab='Survival Status', ylab='Expression Value')
xyplot(values ~ fBatch | ind, data=dfData, type='p', varwidth=F)

##################################
# library(LearnBayes)
logit.inv = function(p) {exp(p)/(exp(p)+1) }

## binomial prediction
mypred = function(theta, data){
  betas = theta # vector of betas i.e. regression coefficients for population
  ## data
  mModMatrix = data$mModMatrix
  # calculate fitted value
  iFitted = mModMatrix %*% betas
  # using logit link so use inverse logit
  iFitted = logit.inv(iFitted)
  return(iFitted)
}
# ## lets write a custom glm using a bayesian approach
# ## write the log posterior function
# mylogpost = function(theta, data){
#   ## parameters to track/estimate
#   betas = theta # vector of betas i.e. regression coefficients for population
#   ## data
#   resp = data$resp # resp
#   mModMatrix = data$mModMatrix
#   
#   # calculate fitted value
#   iFitted = mModMatrix %*% betas
#   # using logit link so use inverse logit
#   iFitted = logit.inv(iFitted)
#   # write the priors and likelihood 
#   lp = dnorm(betas[1], 0, 10, log=T) + sum(dnorm(betas[-1], 0, 10, log=T))
#   lik = sum(dbinom(resp, 1, iFitted, log=T))
#   val = lik + lp
#   return(val)
# }

dfData = dfData.org
dfData = dfData[,colnames(dfData) %in% CVariableSelection.ReduceModel.getMinModel(oVar.sub, 7)]
dim(dfData)
#dfData = scale(dfData)
dim(dfData)
dfData = data.frame(dfData, fGroups)

# lData = list(resp=ifelse(dfData$fGroups == '0', 0, 1), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))
# start = c(rep(0, times=ncol(lData$mModMatrix)))
# mylogpost(start, lData)
# 
# fit.2 = laplace(mylogpost, start, lData)
# fit.2
# data.frame(coef(fit.1), fit.2$mode)
# se = sqrt(diag(fit.2$var))


### lets take a sample from this 
## parameters for the multivariate t density
# tpar = list(m=fit.2$mode, var=fit.2$var*2, df=4)
# ## get a sample directly and using sir (sampling importance resampling with a t proposal density)
# s = sir(mylogpost, tpar, 5000, lData)
# colnames(s) = colnames(lData$mModMatrix)
# apply(s, 2, mean)
# apply(s, 2, sd)
# pairs(s, pch=20)
# fit.2$sir = s

lData = list(resp=ifelse(dfData$fGroups == '0', 0, 1), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))

stanDso = rstan::stan_model(file='binomialRegression.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

## give initial values
initf = function(chain_id = 1) {
  list(betas=rep(0, times=ncol(lStanData$X)), tau=0.5)
}


fit.stan = sampling(stanDso, data=lStanData, iter=2000, chains=2, pars=c('tau', 'betas2'), init=initf, 
                    control=list(adapt_delta=0.99, max_treedepth = 11))

print(fit.stan, c('betas2', 'tau'))
print(fit.stan, 'tau')
traceplot(fit.stan, 'tau')

## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$betas2
dim(mCoef)
colnames(mCoef) = c('Intercept', CVariableSelection.ReduceModel.getMinModel(oVar.sub, 7))
pairs(mCoef, pch=20)

## get the predicted values
dfData.new = dfData
str(dfData.new)
## create model matrix
X = as.matrix(cbind(rep(1, times=nrow(dfData.new)), dfData.new[,colnames(mCoef)[-1]]))
colnames(X) = colnames(mCoef)
ivPredict = mypred(colMeans(mCoef), list(mModMatrix=X))
fPredict = rep('0', times=length(ivPredict))
fPredict[ivPredict > 0.64556962] = '1'
table(fPredict, fGroups)

## fit a binomial model
fit.bin = glm(fGroups ~ ., data=dfData, family='binomial')
summary(fit.bin)
ivPredict.bin = predict(fit.bin, type = 'response')

m = data.frame(round(ivPredict, 2), round(ivPredict.bin, 2), fGroups)

## find the optimal point
ivTruth = fGroups == '1'

p = prediction(ivPredict, ivTruth)
perf.alive = performance(p, 'tpr', 'fpr')
dfPerf.alive = data.frame(c=perf.alive@alpha.values, t=perf.alive@y.values[[1]], f=perf.alive@x.values[[1]], 
                          r=perf.alive@y.values[[1]]/perf.alive@x.values[[1]])
colnames(dfPerf.alive) = c('c', 't', 'f', 'r')


ivTruth = !ivTruth
p = prediction(1-ivPredict, ivTruth)
perf.death = performance(p, 'tpr', 'fpr')
dfPerf.death = data.frame(c=perf.death@alpha.values, t=perf.death@y.values[[1]], f=perf.death@x.values[[1]], 
                          r=perf.death@y.values[[1]]/perf.death@x.values[[1]])
colnames(dfPerf.death) = c('c', 't', 'f', 'r')

plot(perf.alive)
plot(perf.death, add=T, col='red')

dfPlot = data.frame(fGroups, ivPredict)

xyplot(ifelse(fGroups == '1', 1, 0) ~ ivPredict, type=c('p'), ylab='Survived')

plot(ivPredict, ifelse(fGroups == '1', 1, 0), type=c('p'), ylab='Class Prediction', pch=20)
#lines(lowess(dfPlot$ivPredict[dfPlot$tpr > 0.9], ifelse(dfPlot$fGroups[dfPlot$tpr > 0.9] == '1', 1, 0)))
lines(lowess(dfPlot$ivPredict, ifelse(dfPlot$fGroups == '1', 1, 0)))

points(ivPredict, ifelse(fGroups == '0', 1, 0), type=c('p'), col='red', pch=20)
lines(lowess(ivPredict, ifelse(fGroups == '0', 1, 0)), col='red')
plot(1-ivPredict, ifelse(fGroups == '0', 1, 0))

# # ## get the intercept at population level
# iIntercept = mCoef[,1]
# mCoef = mCoef[,-1]
# # ## add the intercept to each coefficient, to get the full coefficient
# # mCoef = sweep(mCoef, 1, iIntercept, '+')
# 
# ## function to calculate statistics for a coefficient
# getDifference = function(ivData){
#   # get the difference vector
#   d = ivData
#   # get the z value
#   z = mean(d)/sd(d)
#   # get 2 sided p-value
#   p = pnorm(-abs(mean(d)/sd(d)))*2
#   return(p)
# }
# 
# ivPval = apply(mCoef, 2, getDifference)
# hist(ivPval)
# plot(colMeans(mCoef), ivPval, pch=19)
# m = colMeans(mCoef)
# names(m) = colnames(lData$mModMatrix)[2:length(m)]
# m = abs(m)
# m = sort(m, decreasing = T)




















