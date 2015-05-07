#Author: Kunal Bhutani
#Date: May 6, 2015
#Purpose: Find changepoints that indicate a shift in metabolite levels

library(changepoint)

# Load data
source("R/metabolomics/load.R")

# Output Directory
local.path = config$paths$local
today <- Sys.Date()
date = format(today, format="%m%d%Y")
analysis.path = file.path(local.path, "analysis", "metabolomics", date)
dir.create(analysis.path)
cp.path = file.path(analysis.path, "change_points")
dir.create(cp.path)

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
                      "samples" = rownames(scaled.metabolite.df),
                      "metabolites" = colnames(not.na.metabolites),
                      "m.pelt"= m.pelt, "m.binseg"= m.binseg,
                      "v.pelt"= v.pelt, "v.binseg"= v.binseg,
                      "mv.pelt"= mv.pelt, "mv.binseg"= mv.binseg)
  changepoints
  
}

PlotChangePoints = function(changepoints, output.prefix, output.suffix, title.prefix,
                            ylabel="Raw Counts") {
  
  PlotMetabolite = function (metabolite) {
    metabolite.id = substr(metabolite, 2, nchar(metabolite))
    biochemical = metab.info[metabolite.id, 'BIOCHEMICAL']
    out.fn = paste(output.prefix, biochemical, output.suffix, sep="")
    jpeg(out.fn)
    labels = sapply(rownames(scaled.metabolite.df),
                    function(x) {tail(strsplit(x, "_")[[1]], n=1)})
    plot(changepoints$m.pelt[[metabolite]], type="b", cpt.col="red", xlab="", ylab="", xaxt="n", yaxt="n")
    par(new=T)
    plot(changepoints$v.pelt[[metabolite]], type="b", cpt.col="blue", xlab="", ylab="", xaxt="n", yaxt="n")
    par(new=T)
    plot(changepoints$v.pelt[[metabolite]], type="b", cpt.col="green", xaxt="n",
         xlab="Date", ylab=ylabel, main=paste(title.prefix, biochemical))
    axis(1, at=1:length(labels), labels=labels, cex.axis=0.95)
    legend("topright", lty=c(1,1,1), lwd=c(2.5,2.5,2.5),col=c("red", "blue", "green"),
           legend=c("Mean Change","Var Changes", "MeanVar Changes")) 
    dev.off()
    
  }
  
  sapply(changepoints$metabolites, PlotMetabolite)
  
}


# Initial analysis only on run 6
run.6 = LoadData(c(6))
metabolite.df = run.6$metabolite.df[,3:dim(run.6$metabolite.df)[2]]
rownames(metabolite.df) = run.6$metabolite.df$Xsample.cols

output.prefix = file.path(cp.path, "run6-")
title.prefix = "Run 6: "
# Scale across samples by log-transforming data
scaled.log.metabolite.df = log(metabolite.df)
changepoints = CalculateChangePoints(scaled.log.metabolite.df)
output.suffix = "-log.jpg"
ylabel = "Log Raw Score"
PlotChangePoints(changepoints, output.prefix, output.suffix, title.prefix, ylabel)

# Using standard scaling with mean = 0, sd = 1 here.
scaled.z.score.metabolite.df = scale(metabolite.df)
changepoints = CalculateChangePoints(scaled.z.score.metabolite.df)
output.suffix = "-zscore.jpg"
ylabel = "Z-Score"
PlotChangePoints(changepoints, output.prefix, output.suffix, title.prefix, ylabel)

# Scale by median
scaled.median.metabolite.df = apply(scaled.metabolite.df, 2, function(x) { x/median(x)})
changepoints = CalculateChangePoints(scaled.median.metabolite.df)
output.suffix = "-median.jpg"
ylabel = "Median Scaled"
PlotChangePoints(changepoints, output.prefix, output.suffix, title.prefix, ylabel)



