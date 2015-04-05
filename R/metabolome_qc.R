#Author: Kunal Bhutani
#Date: March 31st, 2015
#This script does basic QC for the metabolome data. 

# Libraries
library(xlsx)

# Paths
local.path = "/home/kunal/tscc_projects/tanner/data/family3/metabolome"
date = "04042015"
analysis.path = file.path(local.path, "analysis", date)
tn.path = file.path(analysis.path, "technical_noise")
path.to.Rdata = file.path(local.path, "metabolome.data_032715.RData")
rownames(metab.info) = metab.info[,'CHEMICAL.ID']


# Loading in Data
load(path.to.Rdata)

# Functions
GetAllSamples <- function(runs = 1:length(plasma)) {
  Reduce(union, lapply(plasma[runs], function(x) names(x)[grepl("X001", names(x))]))
}


CreateSampleRunsDf <- function(samples) {
  # Currently creates a sample/run df for the plasma samples; can be 
  # modified to be more generalized.
  # Creates a dataframe:
  #  cols: samples
  #  rows: run id
  #  values: column in the run df to which the sample belongs, NA if not in run
  sample.runs.df = data.frame(lapply(plasma, function (x) match(samples, colnames(x))))
  colnames(sample.runs.df) = 1:length(plasma)
  rownames(sample.runs.df) = samples
  sample.runs.df
}


FindOverlappingMetabolites <- function(runs = 1:length(plasma)) {
  # Finds a set of metabolites that are found in each of the runs
  # for the plasma samples
  Reduce(intersect, lapply(plasma[runs], function(x) rownames(x)))
}

CreateSampleDfs <- function(sample.runs.df, metabolites) {
  # Creates a list of sample specific dataframes based on the values in sample.runs.df
  
  ProcessSample <- function(sample) {
    sample.runs = sample.runs.df[sample, ]
    sample.runs = sample.runs[which(!is.na(sample.runs))]
    sample.run.idx = as.numeric(colnames(sample.runs))
    sample.df = data.frame(lapply(sample.run.idx, function(x) plasma[[x]][metabolites,
                                                                          sample]))
    rownames(sample.df) = metabolites
    colnames(sample.df) = sample.run.idx
    sample.df
  }
  sample.dfs = lapply(rownames(sample.runs.df), ProcessSample)
  names(sample.dfs) = rownames(sample.runs.df)
  sample.dfs
}

CreateMetaboliteDfs <- function(sample.runs.df, metabolites) {
  # TODO: Create Metabolite Dfs, makes it easier for later analysis for specific metabolites
}

# Analysis


# TODO:
# Find how many metabolites pass based on correlation test vs coefficient of variation
# Use different samples as the scale
# Create beautiful plots


#sample1 = 'X001_002_Metabolites_Plasma_JCVI.00001_04.16.14'
#sample2 = 'X001_002_Metabolites_Plasma_JCVI.00002_04.16.14'
#samples = c(sample1, sample2)

# Create Sample Specific Dfs
samples = GetAllSamples()
sample.runs.df = CreateSampleRunsDf(samples)
overlapping.metabolites = FindOverlappingMetabolites()
sample.dfs = CreateSampleDfs(sample.runs.df, overlapping.metabolites)

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
rownames(sample1.values) = metab.info[rownames(sample1.values),]$BIOCHEMICAL
rownames(sample2.values) = rownames(sample1.values)
for (i in 1:length(sample1.values)){
  metabolite = rownames(sample1.values)[i]
  fn = file.path(tn.path, sprintf("%s_%s_v_%s", metabolite, sample1, sample2))
  png(fn)
  title = sprintf("%s for JCVI.00001 vs JCVI.00002 \n Corr: %0.2f", metabolite, correlations[metabolite])
  plot(as.numeric(sample1.values[metabolite,]), as.numeric(sample2.values[metabolite,]),
       xlab="Raw Counts JCVI.00001", ylab="Raw Counts JCVI.00002", main=title )
  dev.off()
}

# Calculate Coefficient of Variation based on Scaled for TP2

# Create new scaled values based on sample 1 only.

# Coefficient of Variation: sd/mean

