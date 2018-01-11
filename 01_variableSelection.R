# Name: 01_variableSelection.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 11/01/2018
# Desc: initial report for variable selection using the mass spec data

if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/master/CCrossValidation.R'
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

# calculate a scaling factor
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
sf = rowMeans(mDat)
mDat.res = sweep(mDat, 1, sf, '-')
## use median as size factor
sf = apply(mDat.res, 2, function(x) quantile(x, prob=0.5))
mDat.norm = sweep(mDat, 2, sf, '-')

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
fBatch = fGroups[-iOutliers]

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

## perform analysis with either normalised or raw data
#mDat = mDat.norm
# calculate the coef of var for each gene
cv = apply(mDat, 1, function(x) sd(x)/abs(mean(x)))
# check cv
summary(cv)

# cut data into groups based on quantiles of cv
# cut.pts = quantile(cv, probs = 0:10/10)
# groups = cut(cv, breaks = cut.pts, include.lowest = T)
# groups = cut(cv, breaks = cut.pts, include.lowest = T, labels = 0:9)
# iMean = apply(mDat, 1, mean)
# iVar = apply(mDat, 1, var)
# coplot((cv) ~ iMean | groups)
# coplot(iVar ~ iMean | groups)
# # choose genes with small cv
# # f = cv < 0.2
# # choosing groups from quantile 0 to 40
# mDat.lowcv = mDat[groups %in% c(0),]
# dim(mDat.lowcv)
# 
# sf = rowMeans(mDat.lowcv)
# mDat.res = sweep(mDat.lowcv, 1, sf, '-')
# sf = apply(mDat.res, 2, median)
# mDat.norm2 = sweep(mDat, 2, sf, '-')
# 
# oDiag.1 = CDiagnosticPlots(mDat.norm, 'lneg norm')
# oDiag.2 = CDiagnosticPlots(mDat.norm2, 'lneg norm2')

dim(mDat)
# select a subset of genes that show differential expression
p.t = apply(mDat, 1, function(x) t.test(x ~ fGroups[-iOutliers])$p.value)
p.t.adj = p.adjust(p.t, 'BH')
p.t.adj = sort(p.t.adj, decreasing = F)
t = names(p.t.adj[1:200])
cvTopFeatures.200 = t

mDat = mDat[rownames(mDat) %in% cvTopFeatures.200, ]

## reload the full data
dfData = read.csv('temp/transposed_cleaned.csv', header = T, na.strings = c('na', 'NA', 'NaN'))
dfData = dfData[-iOutliers,-c(1:3)]
dim(dfData)
## add the new normalised/raw data
colnames(dfData)[1:100]
dfData = dfData[,1:90]
dim(mDat)
dim(dfData)
identical(rownames(dfData), colnames(mDat))
dfData = cbind(dfData, t(mDat))
dim(dfData)

gc(reset = T)
## drop covariates with lots of with NAs
s = sapply(1:ncol(dfData), function(x) sum(is.na(dfData[,x])))
i = which(s > 7)
length(i)
s[i]
dfData = dfData[,-i]
dim(dfData)
dfData = na.omit(dfData)
dim(dfData)
fGroups = factor(dfData$X90.day)
dfData = dfData[,-(which(colnames(dfData) %in% c('X90.day', 'X30.day', 'X1.Year', "Sample.File.Name")))]

# save backup
dfData.bk = dfData

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
#cvKeep = c('X1091.89_623.234', 'Neutrophil', 'MELD', 'X1106.76_466.66', 'X1091.88_607.132', 'X939.59_348.571')
cvKeep = names(colSums(s)[colSums(s) <= 3])
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