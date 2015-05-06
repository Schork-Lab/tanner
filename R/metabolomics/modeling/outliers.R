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

# Find outliers

# Method 1: z-score > 2
metabolite.z.scores = scale(scaled.metabolite.df) # scale across metabolites
z.score.outliers = which(abs(scale(scaled.metabolite.df)) > 2, arr.ind=T)[,2]


