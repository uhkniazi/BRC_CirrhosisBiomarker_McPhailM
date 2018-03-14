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
dfData = dfData[,-(which(colnames(dfData) %in% c('Height_', 'Weight', 'Gender_.1.male.', "fGroups",
                                                 'fControl')))]

# save backup
dfData.bk = dfData
fGroups.bk = fGroups
fControl.bk = fControl

############### load the mass spec data lpos
dfData = read.csv('temp/transposed_cleaned_lpos.csv', header = T, na.strings = c('na', 'NA', 'NaN'))

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


