# Name: 05_jointModel.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 14/03/2018
# Desc: combined model for data from clinical and mass spec

source('header.R')

if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/experimental/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

dfData = read.csv('~/Dropbox/Home/temp/transposed_cleaned_lpos.csv', header = T, na.strings = c('na', 'NA', 'NaN'))
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
dfData = dfData[,-(which(colnames(dfData) %in% c('Height_', 'Weight', 'Gender_.1.male.', "fGroups",
                                                 'fControl')))]

# save backup
dfData.bk = dfData
fGroups.bk = fGroups
fControl.bk = fControl

############### load the mass spec data lpos
dfData = read.csv('~/Dropbox/Home/temp/transposed_cleaned_lpos.csv', header = T, na.strings = c('na', 'NA', 'NaN'))

## match the 2 data sets
i = match(dfData.bk$sample_ID_V2, dfData$sample_ID_V2)
dfData = dfData[i,]
dim(dfData)
identical(as.character(dfData$sample_ID_V2), as.character(dfData.bk$sample_ID_V2))

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

## perform analysis with normalised
mDat = mDat.norm
dim(mDat)
rm(mDat.res); rm(mDat.norm)

dfData = dfData.bk
dfData = cbind(dfData, t(mDat))
colnames(dfData)[1:100]
dfData = dfData[,-1]

dfData.org = dfData

##################### choose the covariates based on previous analysis
cvModel.1 = scan(what=character())
cvModel.2 = scan(what=character())

dfData.1 = dfData.org[,cvModel.1]
dfData.2 = dfData.org[,cvModel.2]

############ binomial models for each data set
dfData.1$fGroups = fGroups
fit.1 = glm(fGroups ~ ., data=dfData.1, family = binomial(link='logit'))
summary(fit.1)

dfData.2$fGroups = fGroups
fit.2 = glm(fGroups ~ ., data=dfData.2, family = binomial(link='logit'))
summary(fit.2)

######## mixture model in stan
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='binomialResp2ComponentRegression.stan')

## make 2 model matrices
m1 = model.matrix(fGroups ~ ., data=dfData.1)
m2 = model.matrix(fGroups ~ ., data=dfData.2)
m1 = as.matrix(m1[,-1])
m2 = as.matrix(m2[,-1])
## prepare data
lStanData = list(Ntotal=nrow(dfData.1), Ncol1=ncol(m1), Ncol2=ncol(m2),
                 X1=m1, X2=m2, iMixtures=2, y=ifelse(fGroups == '0', 0, 1))

fit.stan = sampling(stanDso, data=lStanData, iter=50000, chains=1, cores=1)
print(fit.stan, c('betasMix1_2', 'betasMix2_2', 'mu', 'iMixWeights', 'tau'), digi=3)
traceplot(fit.stan, c('betasMix1', 'betasMix2', 'mu', 'iMixWeights'))

########################################## fitted values
## get fitted values
logit.inv = function(p) {exp(p)/(exp(p)+1) }

m = extract(fit.stan)
names(m)
intercepts = apply(m$mu, 2, mean)
betas.1 = colMeans(m$betasMix1_2)
betas.2 = colMeans(m$betasMix2_2)
iMixWeights = colMeans(m$iMixWeights)

m1 = model.matrix(fGroups ~ ., data=dfData.1)
m2 = model.matrix(fGroups ~ ., data=dfData.2)

iPred1 = logit.inv(m1 %*% c(intercepts[1], betas.1))
iPred2 = logit.inv(m2 %*% c(intercepts[2], betas.2))

head(round(cbind(iPred1, iPred2), 2))

library(lattice)
xyplot(iPred1 ~ fGroups, xlab='Actual Group', ylab='Predicted Probability of Being Alive (1)')
xyplot(iPred2 ~ fGroups, xlab='Actual Group', ylab='Predicted Probability of Being Alive (1)')

## get aggregate
iAggregate = cbind(iPred1, iPred2)
iAggregate = sweep(iAggregate, 2, iMixWeights, '*')
iAggregate = rowSums(iAggregate)

xyplot(iAggregate ~ fGroups, xlab='Actual Group', ylab='Predicted Probability of Being Alive (1)')

## choose an appropriate cutoff for accept and reject regions
ivTruth = fGroups == '1'

p = prediction(iAggregate, ivTruth)
perf.alive = performance(p, 'tpr', 'fpr')
dfPerf.alive = data.frame(c=perf.alive@alpha.values, t=perf.alive@y.values[[1]], f=perf.alive@x.values[[1]], 
                          r=perf.alive@y.values[[1]]/perf.alive@x.values[[1]])
colnames(dfPerf.alive) = c('c', 't', 'f', 'r')

ivTruth = !ivTruth
p = prediction(1-iAggregate, ivTruth)
perf.death = performance(p, 'tpr', 'fpr')
dfPerf.death = data.frame(c=perf.death@alpha.values, t=perf.death@y.values[[1]], f=perf.death@x.values[[1]], 
                          r=perf.death@y.values[[1]]/perf.death@x.values[[1]])
colnames(dfPerf.death) = c('c', 't', 'f', 'r')

par(mfrow=c(1,1))
plot(perf.alive)
plot(perf.death, add=T, col='red')
legend('bottomright', legend = c('Alive', 'Dead'), col = 1:2, lty=1)

fPredict = rep('reject', times=length(iAggregate))
fPredict[iAggregate >= 8.329126e-01] = '1'
fPredict[iAggregate <= (1-0.2570353157)] = '0'
table(fPredict, fGroups)

## draw these accept reject points
xyplot(iAggregate ~ fGroups, xlab='Actual Group', ylab='Predicted Probability of Being Alive (1)', groups=fPredict,
       auto.key = list(columns=3))


################## fit one binomial model rather than mixture
##################################
dfData = dfData.org
dfData = dfData[,c(cvModel.1, cvModel.2)]
dim(dfData)
dfData$fGroups = fGroups

lData = list(resp=ifelse(dfData$fGroups == '0', 0, 1), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='binomialRegression2Scales.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

head(lData$mModMatrix)
lStanData$ivCoefMap = c(5, 12)

fit.stan = sampling(stanDso, data=lStanData, iter=2000, chains=2, pars=c('tau', 'betas2', 'betas'), 
                    control=list(adapt_delta=0.99, max_treedepth = 12))

print(fit.stan, c('betas2', 'tau'))
print(fit.stan, 'tau')
traceplot(fit.stan, 'tau')

fit.3 = glm(fGroups ~ ., data=dfData, family = binomial(link='logit'))
summary(fit.3)

## prediction
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


mCoef = extract(fit.stan)$betas2
dim(mCoef)
colnames(mCoef) = c('Intercept', c(cvModel.1, cvModel.2))
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

############################ variable selection section
dfData = dfData.org
dfData = dfData[,c(cvModel.1, cvModel.2)]
dim(dfData)

## create a test and training set
set.seed(1234);
test = sample(1:nrow(dfData), size = nrow(dfData)*0.2, replace = F)
table(fGroups[test])
table(fGroups[-test])


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
for (i in 1:11){
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

dev.off(dev.cur())
