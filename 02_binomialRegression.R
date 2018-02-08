# Name: 02_binomialRegression.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 1/02/2018
# Desc: testing a binomial regression approach

dfData = read.csv(file.choose(), header = T, na.strings = c('na', 'NA', 'NaN'))
dfData = t(dfData)
dfData = data.frame(dfData)
cn = t(dfData[1,])
colnames(dfData) = cn[,1]
dfData = dfData[-1,]
rm(cn)

## remove the unwanted groups 
i = which(dfData$`Class ID` %in% c('ACLF', 'AD', 'SC'))
length(i)
dfData = dfData[i,]
dim(dfData)
dfData = droplevels.data.frame(dfData)

gc(reset = T)

write.csv(dfData, 'temp/transposed_cleaned.csv', row.names = F)
## clears data types etc, do a fresh reload
dfData = read.csv('temp/transposed_cleaned.csv', header = T, na.strings = c('na', 'NA', 'NaN'))
gc(reset = T)
## main grouping factor
fGroups = factor(dfData$X90.day)
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

m = colMeans(mDat)
plot(m)
## 3 samples are different or outliers drop those
k = hclust(dist(m))
plot(k)
c = cutree(k, k = 2)
table(c)
iOutliers = which(c == 2)

## drop the outliers
mDat = mDat[,-iOutliers]
dim(mDat)
dfData = read.csv('temp/transposed_cleaned.csv', header = T, na.strings = c('na', 'NA', 'NaN'))
gc(reset = T)
dfData = dfData[-iOutliers, ]
dim(dfData)
identical(rownames(dfData), colnames(mDat))

gc(reset = T)
fGroups = factor(dfData$X90.day)

dim(mDat)
sf = rowMeans(mDat)
mDat.res = sweep(mDat, 1, sf, '-')
## use median as size factor
sf = apply(mDat.res, 2, function(x) quantile(x, prob=0.5))
mDat = sweep(mDat, 2, sf, '-')

dfData = data.frame(t(mDat), fGroups)

# save backup
dfData.bk = dfData
mDat.bk = mDat

## choose a subset of the data for modelling
dfData = dfData[,cvSample[1]]
mDat = mDat[cvSample[1],]
identical(rownames(dfData), colnames(mDat))
identical(colnames(dfData), rownames(mDat))
dfData = data.frame(t(mDat), fGroups)

### DE model using limma
library(limma)

design = model.matrix(~ dfData$fGroups)
#colnames(design) = levels(lData$cov1)
head(design)

mData = as.matrix(dfData[, -ncol(dfData)])

fit = lmFit(t(mData), design)
fit = eBayes(fit)

dfLimmma.2 = topTable(fit, coef = 2, adjust='BH', number=Inf)

################### random effects model
### utility functions
# utility function to calculate gamma distribution used as a prior for scale
gammaShRaFromModeSD = function( mode , sd ) {
  # function changed a little to return jeffery non-informative prior
  if ( mode <=0 || sd <=0 ) return( list( shape=0.5 , rate=0.0001 ) ) 
  rate = ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
  shape = 1 + mode * rate
  return( list( shape=shape , rate=rate ) )
}

f_plotVolcano = function(dfGenes, main, p.adj.cut = 0.1, fc.lim = c(-3, 3)){
  p.val = -1 * log10(dfGenes$P.Value)
  fc = dfGenes$logFC
  # cutoff for p.value y.axis
  y.cut = -1 * log10(0.01)
  col = rep('lightgrey', times=length(p.val))
  c = which(dfGenes$adj.P.Val < p.adj.cut)
  col[c] = 'red'
  plot(fc, p.val, pch=20, xlab='Fold Change', ylab='-log10 P.Value', col=col, main=main, xlim=fc.lim)
  abline(v = 0, col='grey', lty=2)
  abline(h = y.cut, col='red', lty=2)
  # second cutoff for adjusted p-values
  y.cut = quantile(p.val[c], probs=0.95)
  abline(h = y.cut, col='red')
  # identify these genes
  g = which(p.val > y.cut)
  lab = dfGenes[g, 'SYMBOL']
  text(dfGenes$logFC[g], y = p.val[g], labels = lab, pos=2, cex=0.6)
}

plotMeanFC = function(m, dat, p.cut, title){
  col = rep('grey', length.out=(nrow(dat)))
  col[dat$adj.P.Val < p.cut] = 'red'
  #mh = cut(m, breaks=quantile(m, 0:50/50), include.lowest = T)
  plot(m, dat$logFC, col=col, pch=20, main=title, xlab='log Mean', ylab='logFC', ylim=c(-2, 2), cex=0.6)
}

### end utility functions

dfData = data.frame(t(mDat))
dfData = stack(dfData)
dfData$fBatch = fGroups
dfData$Coef = factor(dfData$fBatch:dfData$ind)
dfData = droplevels.data.frame(dfData)
dfData = dfData[order(dfData$Coef), ]
dim(dfData)

# #### fit mixed effect model
library(lme4)
fit.lme1 = lmer(values ~ 1 + (1 | Coef), data=dfData, REML=F)
summary(fit.lme1)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='normFiniteMixture1RandomEffect.stan')

## calculate hyperparameters for variance of coefficients
l = gammaShRaFromModeSD(sd(dfData$X1000.44_415.827), 2*sd(dfData$X1000.44_415.827))

## set initial values
# ran = ranef(fit.lme1)
# r1 = ran$Coef
# 
# initf = function(chain_id = 1) {
#   list(sigmaRan1 = 2, rGroupsJitter1=r1)
# }

### try a t model without mixture
lStanData = list(Ntotal=nrow(dfData), Nclusters1=nlevels(dfData$fGroups),
                 NgroupMap1=as.numeric(dfData$fGroups),
                 y=dfData$X760.369_26.816, iMixtures=2,
                 gammaShape=l$shape, gammaRate=l$rate, iIntercepts=tapply(dfData$X760.369_26.816, dfData$fGroups, mean))

fit.stan = sampling(stanDso, data=lStanData, iter=500, chains=2,
                    pars=c('sigmaRan1', 'sigma',
                           'rGroupsJitter1', 'mu', 'iMixWeights', 'muFitted'),
                    cores=2, init=initf, control=list(adapt_delta=0.99, max_treedepth = 12))
save(fit.stan, file='temp/fit.stan.normMix.rds')

print(fit.stan, c('sigmaRan1', 'sigma', 'mu', 'iMixWeights'), digits=3)
print(fit.stan, 'rGroupsJitter1')

## check if labelling degeneracy has occured
## see here: http://mc-stan.org/users/documentation/case-studies/identifying_mixture_models.html
params1 = as.data.frame(extract(fit.stan, permuted=FALSE)[,1,])
params2 = as.data.frame(extract(fit.stan, permuted=FALSE)[,2,])

## check if the means from different chains overlap
## Labeling Degeneracy by Enforcing an Ordering
par(mfrow=c(1,2))
plot(params1$`mu[1]`, params1$`mu[2]`, pch=20, col=2)
plot(params2$`mu[1]`, params2$`mu[2]`, pch=20, col=3)

par(mfrow=c(1,1))
plot(params1$`mu[1]`, params1$`mu[2]`, pch=20, col=2)
points(params2$`mu[1]`, params2$`mu[2]`, pch=20, col=3)

## get the coefficient of interest - Modules in our case from the random coefficients section
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

### compare the results from the 3 models
# dfResults = dfResults[order(dfResults$pvalue, decreasing = F), ]
# dfLimmma.2 = dfLimmma.2[order(dfLimmma.2$P.Value, decreasing = F), ]
dfResults$logFC = dfResults$difference
dfResults$P.Value = dfResults$pvalue
dfLimmma.2$SYMBOL = as.character(rownames(dfLimmma.2))
dfResults$SYMBOL = as.character(rownames(dfResults))

## produce the plots 
f_plotVolcano(dfLimmma.2, 'limma', fc.lim = c(-2, 2))
f_plotVolcano(dfResults, 'Stan', fc.lim=c(-2.5, 2.5))

m = tapply(dfData$values, dfData$ind, mean)
i = match(rownames(dfResults), names(m))
m = m[i]
identical(names(m), rownames(dfResults))
plotMeanFC(m, dfResults, 0.1, 'Stan 2M vs 12M')

m = tapply(dfData$values, dfData$ind, mean)
i = match(rownames(dfLimmma.2), names(m))
m = m[i]
m = m[!is.na(m)]
i = match(names(m), rownames(dfLimmma.2))
dfLimmma.2 = dfLimmma.2[i,]
identical(names(m), rownames(dfLimmma.2))
plotMeanFC(m, dfLimmma.2, 0.1, 'limma 2M vs 12M')


i = match(rownames(dfResults), rownames(dfLimmma.2))
dfLimmma.2 = dfLimmma.2[i,]
i = match(rownames(dfLimmma.2), rownames(dfResults))
dfResults = dfResults[i,]
identical(rownames(dfResults), rownames(dfLimmma.2))

plot(dfResults$pvalue, dfLimmma.2$P.Value, pch=20, cex=0.6, col='grey', main='P Values 2M vs 12M', xlab='Stan', ylab='Limma')
abline(lm(dfLimmma.2$P.Value ~ dfResults$pvalue), col=2, lwd=2)
plot(dfResults$logFC, dfLimmma.2$logFC, pch=20, cex=0.8, col='grey', main='Log FC 2M vs 12M', xlab='Stan', ylab='Limma')
abline(lm(dfLimmma.2$logFC ~ dfResults$logFC), col=2, lwd=1)
df = cbind(stan=dfResults$pvalue, limma=dfLimmma.2$P.Value)

write.csv(dfResults, file='temp/stan_t.csv', row.names = F)
write.csv(dfLimmma.2, file='temp/limma.csv', row.names = F)















p.t = apply(dfData, 2, function(x) t.test(x ~ fGroups)$p.value)
p.t.adj = p.adjust(p.t, 'BH')
p.t.adj = sort(p.t.adj, decreasing = F)
t = names(p.t.adj[1:2000])
cvTopFeatures.2000 = t

dfData = dfData[,cvTopFeatures.2000]
dfData = scale(dfData)
dim(dfData)
dfData = data.frame(dfData, fGroups)
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


fit.stan = sampling(stanDso, data=lStanData, iter=500, chains=2, pars=c('tau', 'betas2'), init=initf, 
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
i = which(ivPval < 0.7)
colnames(lData$mModMatrix)[i+1]
colMeans(mCoef)[i]

dfPval = data.frame(do.call(rbind, lPval))






## create a test and training set
test = sample(1:nrow(dfData), size = nrow(dfData)*0.2, replace = F)
table(fGroups[test])
## perform nested random forest
## adjust boot.num as desired
oVar.r = CVariableSelection.RandomForest(dfData[-test, ], fGroups[-test], boot.num = 100, big.warn = F)
# plot the top 20 variables based on importance scort with 95% confidence interval for standard error
par(mfrow=c(1,1))
plot.var.selection(oVar.r)
# get the variables
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
# select the top 30 variables
cvTopGenes = rownames(dfRF)[1:30]

# use the top 30 features to find top combinations of genes
dfData = dfData[,colnames(dfData) %in% cvTopGenes]

## look at colinear variables
m = NULL;

for (i in 1:ncol(dfData)){
  m = cbind(m, dfData[-test ,i])
}
colnames(m) = colnames(dfData)
mCor = cor(m, use="na.or.complete")
library(caret)
### find the columns that are correlated and should be removed
n = findCorrelation((mCor), cutoff = 0.7, names=T)
data.frame(n)
sapply(n, function(x) {
  (abs(mCor[,x]) >= 0.7)
})
s = sapply(n, function(x) {
  (abs(mCor[,x]) >= 0.7)
})
colSums(s)
cvKeep = names(colSums(s)[colSums(s) <= 3])
#cvKeep = c('SOFA', 'MELD', 'Albumin', 'BMI', 'Neutrophil')
n = n[!(n%in% cvKeep)]
i = which(colnames(dfData) %in% n)
cn = colnames(dfData)[-i]

dfData.bk2 = dfData
dfData = dfData[,cn]

oVar.sub = CVariableSelection.ReduceModel(dfData[-test, ], fGroups[-test], boot.num = 100)
# plot the number of variables vs average error rate
plot.var.selection(oVar.sub)

# print variable combinations
for (i in 1:6){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  cat('Variable Count', i, paste(cvTopGenes.sub), '\n')
  #print(cvTopGenes.sub)
}

## 10 fold nested cross validation with various variable combinations
par(mfrow=c(2,2))
# try models of various sizes with CV
for (i in 1:6){
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
##################################