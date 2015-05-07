#Author: Kunal Bhutani
#Date: May 6, 2015
#Purpose: Find changepoints that indicate a shift in metabolite levels

library(changepoint)

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
# scaled.metabolite.df = t(apply(run.6.metabolite.df, 1, function(x) { scale(log(x), scale=F)}))
# colnames(scaled.metabolite.df) = colnames(run.6.metabolite.df)

# Since these samples are run on the same run, no need to scale the data. A simple log transformation
# is sufficient.
scaled.metabolite.df = log(run.6.metabolite.df)

# Change point analysis
# http://www.jstatsoft.org/v58/i03/paper

# Based on Means
not.na.metabolites = scaled.metabolite.df[,which(colSums(is.na(scaled.metabolite.df)) == 0)]
m.pelt = cpt.mean(t(not.na.metabolites), method="PELT")
m.binseg = cpt.mean(t(not.na.metabolites), method="BinSeg")
m.changepoints = data.frame(mean.PELT = I(sapply(m.pelt, cpts)),
                            mean.BinSeg = I(sapply(m.binseg, cpts)))

# Based on Variance
v.pelt = cpt.var(t(not.na.metabolites), method="PELT")
v.binseg = cpt.var(t(not.na.metabolites), method="BinSeg")
v.changepoints = data.frame(var.PELT = I(sapply(v.pelt, cpts)),
                            var.BinSeg = I(sapply(v.binseg, cpts)))

# Based on MeanVar
mv.pelt = cpt.meanvar(t(not.na.metabolites), method="PELT")
mv.binseg = cpt.meanvar(t(not.na.metabolites), method="BinSeg")
mv.changepoints = data.frame(meanvar.PELT = I(sapply(mv.pelt, cpts)),
                             meanvar.BinSeg = I(sapply(mv.binseg, cpts)))

changepoints = cbind(m.changepoints, v.changepoints, mv.changepoints)
