#Author: Kunal Bhutani
#Date: May 6, 2015
#Purpose: Detect metabolites that are deviating from baseline for an individual

# Packages
library(pracma)
library(extremevalues)

# Load data
source("R/metabolomics/load.R")

CalculateUnivariateOutliers = function(scaled.metabolite.df){
  # Calculates Univariate Outliers based on several methods
  # Method 1: Assumes metabolites are normally distributed, and finds samples where abs(z-score) > 2
  # Method 2: Uses extremevalues package to find samples that deviate from expected distribution
  # Method 3: Uses Median Absolute Deviation to find samples that are higher or lower than expected.
  # 
  # Args: scaled.metabolite.df -- a scaled version of the metabolite values.
  #       Note: if multiple runs are included, make sure to scale by runs; otherwise a log transform
  #             might be okay.
  # Returns: outliers dataframe with outliers as detected by the 3 different methods
  
  
  # Method 1: z-score > 2
  metabolite.z.scores = scale(scaled.metabolite.df) # scale across metabolites
  z.score.outliers = data.frame(lower.z.scores = I(apply(metabolite.z.scores, 2, function (x) {which(x < -2)})),
                                higher.z.scores = I(apply(metabolite.z.scores, 2, function (x) {which(x > 2)})))
  
  # Method 2: Extreme values 
  # http://www.cbs.nl/NR/rdonlyres/21A8D00F-E20B-43B2-A95D-3D089833EED3/0/201003x10pub.pdf
  # TODO: Unclear how the different scaling affects these calculations. Ignoring for now.
  not.na.metabolites = abs(scaled.metabolite.df[,which(colSums(is.na(scaled.metabolite.df)) == 0)])
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
  
  outliers
}


# Initial analysis only on run 6
run.6 = LoadData(c(6))
scaled.metabolite.df = run.6$metabolite.df[,3:dim(run.6$metabolite.df)[2]]
rownames(scaled.metabolite.df) = run.6$metabolite.df$Xsample.cols

# Scale across samples by log-transforming data and then subtracting out mean
#scaled.metabolite.df = t(apply(run.6$metabolite.df, 1, function(x) { scale(log(x), scale=F)}))

# Since these samples are run on the same run, no need to scale the data. A simple log transformation
# is sufficient.
scaled.metabolite.df = log(scaled.metabolite.df)
outliers = CalculateUnivariateOutliers(scaled.metabolite.df)

