#Author: Kunal Bhutani
#Date: March 31st, 2015
#This script does basic QC for the metabolome data. 

# Libraries
library(xlsx)

# Load data
source("load.R")

# Set paths
local.path = config$paths$local
code.path = file.path(local.path, 'code/R/')
setwd(code.path)

# Output paths
today <- Sys.Date()
date = format(today, format="%m%d%Y")
analysis.path = file.path(local.path, "analysis", "metabolomics", date)
dir.create(analysis.path)
tn.path = file.path(analysis.path, "technical_noise")
dir.create(tn.path)
corr.path = file.path(tn.path, "correlation")
dir.create(corr.path)
cv.path = file.path(tn.path, "cv")
dir.create(cv.path)

# Find general statistics about runs and individuals
individuals = sapply(rownames(sample.runs.df), function(x) substr(x, 6, 8))
run.info = !is.na(sample.runs.df[names(individuals[order(individuals)]),])
general.fn = file.path(analysis.path, "general.xlsx")
individual.info = table(individuals)
write.xlsx(individual.info, general.fn, sheetName = "Individual Info")
write.xlsx(run.info, general.fn, sheetName = "Run Info",
           append = TRUE)


# Find the correlation for each metabolite between sample 1 and 2
FindCorrelation <- function(sample1.df, sample2.df) {
  overlapping.runs = intersect(colnames(sample1.df),
                               colnames(sample2.df))
  correlation = diag(cor(t(as.matrix(sample1.df[,overlapping.runs])),
                         t(as.matrix(sample2.df[,overlapping.runs]))))
  correlation
}

sample1 = 'X001_002_Metabolites_Plasma_JCVI.00001_04.16.14'
sample2 = 'X001_002_Metabolites_Plasma_JCVI.00002_04.16.14'
correlations = FindCorrelation(sample.dfs[[sample1]], sample.dfs[[sample2]])
names(correlations) = metab.info[names(correlations),]$BIOCHEMICAL

corr.fn = file.path(tn.path, sprintf("corr_%s_v_%s", sample1, sample2))
png(corr.fn)
title = sprintf("JCVI.00001 vs JCVI.00002 \n Correlation for Metabolites \n (%d)", 
                length(correlations))
hist(correlations, main=title, cex=0.5)
dev.off()

# Create metabolite specific scatter plots for ones that have low correlation
# TODO: Clean up code, use metabolite specific dfs.
sample1.df = sample.dfs[[sample1]]
sample2.df = sample.dfs[[sample2]]
overlapping.runs = intersect(colnames(sample1.df),
                             colnames(sample2.df))
 
sample1.values = sample1.df[overlapping.metabolites[which(correlations < 0.5)],
                                       overlapping.runs]
sample2.values = sample2.df[overlapping.metabolites[which(correlations < 0.5)],
                                       overlapping.runs]
rownames(sample1.values) = metab.info[rownames(sample1.values),'BIOCHEMICAL']
rownames(sample2.values) = rownames(sample1.values)
for (i in 1:length(sample1.values)){
  metabolite = rownames(sample1.values)[i]
  fn = file.path(corr.path, sprintf("%s---%s---v---%s.png", metabolite, sample1, sample2))
  png(fn)
  title = sprintf("%s for JCVI.00001 vs JCVI.00002 \n Corr: %0.2f", metabolite, correlations[metabolite])
  plot(as.numeric(sample1.values[metabolite,]), as.numeric(sample2.values[metabolite,]),
       xlab="Raw Counts JCVI.00001", ylab="Raw Counts JCVI.00002", main=title )
  dev.off()
}

# Create new scaled values based on sample 1 only.
ScaleDfs <- function(sample.df, scale.by.df) {
  overlapping.runs = intersect(colnames(sample.df),
                               colnames(scale.by.df))
  sample.df = sample.df[,overlapping.runs]/scale.by.df[,overlapping.runs]
  sample.df
}
scaled.dfs = sapply(sample.dfs, FUN=ScaleDfs, scale.by.df=sample.dfs[[sample1]])

# Calculate coefficient of Variation: sd/mean
cv.dfs = lapply(scaled.dfs, function(x) {apply(as.matrix(x), 1, 
                                               function(y) {sd(y)/mean(y)})})

# Find metabolites with high c.v. for TP2
high.cv.metabolites = which(cv.dfs[[sample2]] > 0.5)

# Find metabolites with over a fold change between TP1 and TP2 and high c.v.
fold.change.metabolites = which(apply(as.matrix(scaled.dfs[[sample2]]), 1,
                                      function (x) {max(x) > 2}))

# Low confidence metabolites are defined as those which have a high c.v. and a high fold change in at
# least one run
low.confidence.metabolites = intersect(names(high.cv.metabolites), names(fold.change.metabolites))

# Using for loop to create figures again
correlations = FindCorrelation(sample.dfs[[sample1]], sample.dfs[[sample2]])
for (i in 1:length(low.confidence.metabolites)) {
  metabolite.id = low.confidence.metabolites[i]
  metabolite = metab.info[metabolite.id,'BIOCHEMICAL']
  fn = file.path(cv.path, sprintf("%s---%s---scaledby---%s.png", metabolite, sample1, sample2))
  png(fn)
  par(mfrow=c(1,2))
  title = sprintf("Scaled JCVI.00002 \n %s \n CV: %0.1f, Corr: %0.2f",
                  metabolite, cv.dfs[[sample2]][metabolite.id], correlations[metabolite.id])
  plot(overlapping.runs, as.numeric(scaled.dfs[[sample2]][metabolite.id, ]),
       main=title, xlab="Run Id", ylab="Scaled against JCVI.00001")
  plot(as.numeric(sample.dfs[[sample1]][metabolite.id, overlapping.runs]),
       as.numeric(sample.dfs[[sample2]][metabolite.id, overlapping.runs]),
       xlab="Raw JCVI.00001", ylab="Raw JCVI.00002", main=metabolite)
  dev.off()
}


# Output Correlation and CVs to Excel Sheets
sample2.corr = FindCorrelation(sample.dfs[[sample1]], sample.dfs[[sample2]])
sample2.cvs = cv.dfs[[sample2]]
sample2.out = data.frame(Correlation=sample2.corr, CV=sample2.cvs)
rownames(sample2.out) = metab.info[rownames(sample2.out), 'BIOCHEMICAL']
corr.fn = file.path(analysis.path, "correlation.and.cv.xlsx")
write.xlsx(sample2.out, corr.fn, sheetName = "Correlation and C.V. for Sample 1 v Sample2")


