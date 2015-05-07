#Author: Kunal Bhutani
#Date: May 6, 2015
#Purpose: Find changepoints that indicate a shift in metabolite levels

library(changepoint)

# Load data
source("R/metabolomics/load.R")

CalculateChangePoints = function(scaled.metabolite.df) {
  # Change point analysis
  # http://www.jstatsoft.org/v58/i03/paper
  
  # Based on Means
  not.na.metabolites = scaled.metabolite.df[,which(colSums(is.na(scaled.metabolite.df)) == 0)]
  m.pelt = cpt.mean(t(not.na.metabolites), method="PELT")
  names(m.pelt) = colnames(not.na.metabolites)
  m.binseg = cpt.mean(t(not.na.metabolites), method="BinSeg")
  names(m.binseg) = colnames(not.na.metabolites)
  m.changepoints = data.frame(mean.PELT = I(sapply(m.pelt, cpts)),
                              mean.BinSeg = I(sapply(m.binseg, cpts)))
  
  # Based on Variance
  v.pelt = cpt.var(t(not.na.metabolites), method="PELT")
  names(v.pelt) = colnames(not.na.metabolites)
  v.binseg = cpt.var(t(not.na.metabolites), method="BinSeg")
  names(v.binseg) = colnames(not.na.metabolites)
  v.changepoints = data.frame(var.PELT = I(sapply(v.pelt, cpts)),
                              var.BinSeg = I(sapply(v.binseg, cpts)))
  
  # Based on MeanVar
  mv.pelt = cpt.meanvar(t(not.na.metabolites), method="PELT")
  names(mv.pelt) = colnames(not.na.metabolites)
  mv.binseg = cpt.meanvar(t(not.na.metabolites), method="BinSeg")
  names(mv.binseg) = colnames(not.na.metabolites)
  mv.changepoints = data.frame(meanvar.PELT = I(sapply(mv.pelt, cpts)),
                               meanvar.BinSeg = I(sapply(mv.binseg, cpts)))
  
  changepoints = cbind(m.changepoints, v.changepoints, mv.changepoints)
  changepoints = list("changepoints"= changepoints,
                      "m.pelt"= m.pelt, "m.binseg"= m.binseg,
                      "v.pelt"= v.pelt, "v.binseg"= v.binseg,
                      "mv.pelt"= mv.pelt, "mv.binseg"= mv.binseg)
  changepoints
  
}

# Initial analysis only on run 6
run.6 = LoadData(c(6))
scaled.metabolite.df = run.6$metabolite.df[,3:dim(run.6$metabolite.df)[2]]
rownames(scaled.metabolite.df) = run.6$metabolite.df$Xsample.cols

# Scale across samples by log-transforming data and then subtracting out mean
#scaled.metabolite.df = t(apply(scaled.metabolite.df, 1, function(x) { scale(log(x), scale=F)}))

# Since these samples are run on the same run, no need to scale the data. A simple log transformation
# is sufficient.
scaled.metabolite.df = log(scaled.metabolite.df)
changepoints = CalculateChangePoints(scaled.metabolite.df)
