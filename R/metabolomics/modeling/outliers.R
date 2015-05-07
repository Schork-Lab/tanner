#Author: Kunal Bhutani
#Date: May 6, 2015
#Purpose: Detect metabolites that are deviating from baseline for an individual

library(pracma)
library(extremevalues)

# Load data
source("R/metabolomics/load.R")


# Set paths
local.path = config$paths$local
code.path = file.path(local.path, 'code/R/')

# Extract only Run 6
samples = GetAllSamples()
sample.runs.df = CreateSampleRunsDf(samples)
run.6.sample.runs.df = sample.runs.df[,"6",drop=FALSE]
run.6.sample.runs.df = run.6.sample.runs.df[which(!is.na(run.6.sample.runs.df)), , drop=F]
run.6.metabolites = FindOverlappingMetabolites(c(6))
run.6.sample.dfs = CreateSampleDfs(run.6.sample.runs.df, run.6.metabolites)
run.6.metabolite.df = CreateMetaboliteDfs(run.6.sample.runs.df, run.6.sample.dfs,run.6.metabolites)
rownames(run.6.metabolite.df) = run.6.metabolite.df$Xsample.cols
run.6.metabolite.df = run.6.metabolite.df[,3:dim(run.6.metabolite.df)[2]]

# Scale across samples by log-transforming data and then subtracting out mean
#scaled.metabolite.df = t(apply(run.6.metabolite.df, 1, function(x) { scale(log(x), scale=F)}))

# Since these samples are run on the same run, no need to scale the data. A simple log transformation
# is sufficient.
scaled.metabolite.df = log(run.6.metabolite.df)

# Univariate Outliers

# Method 1: z-score > 2
metabolite.z.scores = scale(scaled.metabolite.df) # scale across metabolites
z.score.outliers = data.frame(lower.z.scores = I(apply(metabolite.z.scores, 2, function (x) {which(x < -2)})),
                              higher.z.scores = I(apply(metabolite.z.scores, 2, function (x) {which(x > 2)})))

# Method 2: Extreme values 
# http://www.cbs.nl/NR/rdonlyres/21A8D00F-E20B-43B2-A95D-3D089833EED3/0/201003x10pub.pdf
# TODO: Unclear how the different scaling affects these calculations. Ignoring for now.
not.na.metabolites = abs(scaled.metabolite.df[,which(colSums(is.na(scaled.metabolite.df)) < 3)])
extreme.metabolites = sapply(colnames(not.na.metabolites),  
                                  function(x) { 
                                    outliers = getOutliers(not.na.metabolites[!is.na(not.na.metabolites[,x])
                                                                              ,x],
                                                           method="I")
                                    outliers
                                  })
extreme.outliers = as.data.frame(extreme.metabolites)[c('R2','iLeft','iRight'),]
rownames(extreme.outliers) = c('R2-QQ', 'lower.extreme.values', 'higher.extreme.values')
extreme.outliers = t(extreme.outliers)

# Method 3: Median Absolute Deviation
# http://www.r-bloggers.com/absolute-deviation-around-the-median/
lower.mad.outliers = apply(scaled.metabolite.df, 2, 
                         function(x) {
                           median = median(x, na.rm=T)
                           mad = mad(x, na.rm=T)
                           lower.bound.mad = median + (-2*mad)
                           which(x < lower.bound.mad)
                         })

higher.mad.outliers = apply(scaled.metabolite.df, 2, 
                         function(x) {
                           median = median(x, na.rm=T)
                           mad = mad(x, na.rm=T)
                           upper.bound.mad = median + (2*mad)
                           which(x > upper.bound.mad)
                         })

mad.outliers = data.frame(lower.mad = I(lower.mad.outliers),
                          higher.mad = I(higher.mad.outliers))


# TODO: Method 4: ARIMA Models
# http://stackoverflow.com/questions/13327373/univariate-outlier-detectionex
# http://faculty.chicagobooth.edu/ruey.tsay/teaching/mtsbk/
# http://www.jalobe.com:8080/doc/tsoutliers.pdf
# Although we do not have enough points to get reasonable estimates from these models.


# Combine the outliers
outliers = cbind(z.score.outliers, mad.outliers)
outliers = merge(x=outliers, y=extreme.outliers,
                 by="row.names", all=T)
rownames(outliers) = outliers$Row.names
outliers = outliers[,2:dim(outliers)[2]]