#Author: Kunal Bhutani
#Date: May 6, 2015
#Purpose: Detect metabolites that are deviating from baseline for an individual

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
scaled.metabolite.df = t(apply(run.6.metabolite.df, 1, function(x) { scale(log(x), scale=F)}))
colnames(scaled.metabolite.df) = colnames(run.6.metabolite.df)

# Univariate Outliers

# Method 1: z-score > 2
metabolite.z.scores = scale(scaled.metabolite.df) # scale across metabolites
lower.z.score.outliers = which(metabolite.z.scores < -2, arr.ind=T)
higher.z.score.outliers = which(metabolite.z.scores > 2, arr.ind=T)

# Method 2: Extreme values 
# TODO: Unclear how the different scaling affects these calculations. Ignoring for now.
not.na.metabolites = abs(scaled.metabolite.df[,which(colSums(is.na(scaled.metabolite.df)) < 3)])
extreme.metabolites = sapply(colnames(not.na.metabolites),  
                                  function(x) { 
                                    outliers = getOutliers(not.na.metabolites[!is.na(not.na.metabolites[,x])
                                                                              ,x],
                                                           method="I")
                                    outliers
                                  })

# Method 3: Hampel Filter
# TODO: Error prone, because of the selection of the window
hampel.metabolites = sapply(colnames(not.na.metabolites),  
                             function(x) { 
                               metabolite.levels = not.na.metabolites[!is.na(not.na.metabolites[,x]),x]
                               outliers = hampel(metabolite.levels, 1, 0)
                               outliers
                             })

# Method 4: Median Absolute Deviation
metabolite.mad = apply(scaled.metabolite.df, 2, mad, na.rm=T)
metabolite.median = apply(scaled.metabolite.df, 2, median, na.rm=T)
lower.bound.mad = metabolite.median + (-2*metabolite.mad)
upper.bound.mad = metabolite.median + (2*metabolite.mad)
metabolite.lower = apply(scaled.metabolite.df, 2, 
                         function(x) {
                           median = median(x, na.rm=T)
                           mad = mad(x, na.rm=T)
                           lower.bound.mad = median + (-2*mad)
                           x < lower.bound.mad
                         })

metabolite.higher = apply(scaled.metabolite.df, 2, 
                         function(x) {
                           median = median(x, na.rm=T)
                           mad = mad(x, na.rm=T)
                           upper.bound.mad = median + (2*mad)
                           x > upper.bound.mad
                         })

lower.mad.outliers = which(metabolite.lower, arr.ind=T)
higher.mad.outliers = which(metabolite.higher, arr.ind=T)

