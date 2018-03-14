# Name: 04_variableSelectionClinical.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 14/03/2018
# Desc: variable selection using the clinical data

source('header.R')

if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/experimental/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

dfData = read.csv('temp/transposed_cleaned_lpos.csv', header = T, na.strings = c('na', 'NA', 'NaN'))
gc(reset = T)
## main grouping factor
fGroups = factor(dfData$X90_day)
levels(fGroups)

colnames(dfData)[1:100]

## the first 94 columns are the clinical data
dfData = dfData[,1:94]
## remove garbage
dfData = dfData[,-c(1,2)]
dim(dfData)

gc(reset = T)
## control variable
fControl = factor(dfData$Class_ID)
dfData = dfData[,-(which(colnames(dfData) %in% c('Class_ID', 'X90_day', 'X30_day', 'X1_Year', "Sample.File.Name")))]

# save backup
dfData.bk = dfData
fGroups.bk = fGroups
fControl.bk = fControl

dfData$fGroups = fGroups
dfData$fControl = fControl

## drop covariates with lots of with NAs
s = sapply(1:ncol(dfData), function(x) sum(is.na(dfData[,x])))
i = which(s > 3)
length(i)
s[i]
dfData = dfData[,-i]
dim(dfData)
dfData = na.omit(dfData)
dim(dfData)
dfData = droplevels.data.frame(dfData)
fGroups = dfData$fGroups
fControl = dfData$fControl
str(dfData)
dfData = dfData[,-(which(colnames(dfData) %in% c('sample_ID_V2', 'Height_', 'Weight', 'Gender_.1.male.', "fGroups",
                                                 'fControl')))]

# save backup
dfData.bk = dfData
fGroups.bk = fGroups
fControl.bk = fControl

## create a test and training set
set.seed(1234);
test = sample(1:nrow(dfData), size = nrow(dfData)*0.2, replace = F)
table(fGroups[test])
table(fGroups[-test])

table(fControl[test])
table(fControl[-test])

cvTopVariables = colnames(dfData)

# ## perform nested random forest
# ## adjust boot.num as desired
# oVar.r = CVariableSelection.RandomForest(dfData[-test, cvTopVariables], fGroups[-test], boot.num = 100, big.warn = F)
# # plot the top 20 variables based on importance scort with 95% confidence interval for standard error
# par(mfrow=c(1,1))
# plot.var.selection(oVar.r)
# # get the variables
# dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
# # select the top 30 variables
# cvTopGenes = rownames(dfRF)[1:30]
# 
# # use the top 30 features to find top combinations of genes
# dfData = dfData[,colnames(dfData) %in% cvTopGenes]
# dim(dfData)
# 
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
# cvKeep = names(colSums(s)[colSums(s) <= 2])
# cvKeep = c(cvKeep, 'MELD')
# n = n[!(n%in% cvKeep)]
# i = which(colnames(dfData) %in% n)
# cn = colnames(dfData)[-i]
# 
# dfData.bk2 = dfData
# dfData = dfData[,cn]
# dim(dfData)
# 
# oVar.sub = CVariableSelection.ReduceModel(dfData[-test, ], fGroups[-test], boot.num = 30)
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
dfData.org = dfData.bk
dfData = dfData.bk
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


fit.stan = sampling(stanDso, data=lStanData, iter=2000, chains=2, pars=c('tau', 'betas2'), init=initf, 
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
names(m) = colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)]
m = abs(m)
m = sort(m, decreasing = T)
barplot(m, xlab='', ylab='', main='Ordered Coefficients', xaxt='n')

cvTopGenes.binomial = names(m[m > 0.25])

################### repeat the selection process again 
cvTopGenes = cvTopGenes.binomial;
length(cvTopGenes)
# use the top features to find top combinations of genes
dfData = dfData.org
dfData = dfData[,colnames(dfData) %in% cvTopGenes]
dim(dfData)

oVar.sub = CVariableSelection.ReduceModel(dfData[-test, ], fGroups[-test], boot.num = 100)
# plot the number of variables vs average error rate
plot.var.selection(oVar.sub)

# print variable combinations
for (i in 1:8){
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
dfData = dfData[,colnames(dfData) %in% CVariableSelection.ReduceModel.getMinModel(oVar.sub, 8)]
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
bwplot(log(values) ~ fBatch | ind, data=dfData, type='b', panel=panel.violin, varwidth=F, xlab='Survival Status', ylab='Expression Value')
xyplot(log(values) ~ fBatch | ind, data=dfData, type='p', varwidth=F)

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

dfData = dfData.org
dfData = dfData[,colnames(dfData) %in% CVariableSelection.ReduceModel.getMinModel(oVar.sub, 4)]
dim(dfData)
#dfData = scale(dfData)
dim(dfData)
dfData = data.frame(dfData, fGroups)

lData = list(resp=ifelse(dfData$fGroups == '0', 0, 1), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))

stanDso = rstan::stan_model(file='binomialRegression.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

## give initial values
initf = function(chain_id = 1) {
  list(betas=rep(0, times=ncol(lStanData$X)), tau=0.5)
}


fit.stan = sampling(stanDso, data=lStanData, iter=3000, chains=3, pars=c('tau', 'betas2'), init=initf, 
                    control=list(adapt_delta=0.99, max_treedepth = 11))

print(fit.stan, c('betas2', 'tau'))
print(fit.stan, 'tau')
traceplot(fit.stan, 'tau')

## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$betas2
dim(mCoef)
colnames(mCoef) = c('Intercept', CVariableSelection.ReduceModel.getMinModel(oVar.sub, 4))
pairs(mCoef, pch=20)

## get the predicted values
dfData.new = dfData
str(dfData.new)
## create model matrix
X = as.matrix(cbind(rep(1, times=nrow(dfData.new)), dfData.new[,colnames(mCoef)[-1]]))
colnames(X) = colnames(mCoef)
ivPredict = mypred(colMeans(mCoef), list(mModMatrix=X))
xyplot(ivPredict ~ fGroups, xlab='Actual Group', ylab='Predicted Probability of Being Alive (1)')
## choose an appropriate cutoff for accept and reject regions
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

par(mfrow=c(1,1))
plot(perf.alive)
plot(perf.death, add=T, col='red')
legend('bottomright', legend = c('Alive', 'Dead'), col = 1:2, lty=1)

fPredict = rep('reject', times=length(ivPredict))
fPredict[ivPredict >= 0.958387352] = '1'
fPredict[ivPredict <= (1-7.015914e-01)] = '0'
table(fPredict, fGroups)

## draw these accept reject points
xyplot(ivPredict ~ fGroups, xlab='Actual Group', ylab='Predicted Probability of Being Alive (1)', groups=fPredict,
       auto.key = list(columns=3))

## fit a binomial model
fit.bin = glm(fGroups ~ ., data=dfData, family='binomial')
summary(fit.bin)
ivPredict.bin = predict(fit.bin, type = 'response')