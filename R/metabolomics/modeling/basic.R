#Author: Kunal Bhutani
#Date: April 27th, 2015
#This script does basic modeling for the different metabolites.

# Libraries
library(xlsx)

# Load data
code.path = "/home/kunal/tscc_projects/tanner/code/tanner_project/R/"
setwd(code.path)
source("load.R")

# Output paths
today <- Sys.Date()
date = format(today, format="%m%d%Y")
analysis.path = file.path(local.path, "analysis", date)
dir.create(analysis.path)
spaghetti.plots.path = file.path(analysis.path, "spaghetti")
dir.create(spaghetti.plots.path)


# Scale by a metabolite
metabolite.df[,3:dim(metabolite.df)[2]] = metabolite.df[,3:dim(metabolite.df)[2]]/metabolite.df[,3]

# Output spaghetti plots
for (metabolite in overlapping.metabolites) {
  metabolite.col.name = paste("X", metabolite, sep="")
  chemical.id = metab.info[metabolite,'BIOCHEMICAL']
  fn = file.path(spaghetti.plots.path, paste(chemical.id,'.png',sep=""))
  png(fn)
  interaction.plot(metabolite.df$Xrun.cols,metabolite.df$Xsample.cols,
                   metabolite.df[, metabolite.col.name],, 
                   xlab="Run", ylab="Level", col=c(1:length(sample.dfs)), legend=F,
                   main=chemical.id)
  dev.off()
}

