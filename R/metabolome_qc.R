#Author: Kunal Bhutani
#Date: March 31s, 2015
#This script does basic QC for the metabolome data. 

# Libraries
library(xlsx)

# Paths
local.path = "/home/kunal/tscc_projects/tanner/data/family3/metabolome"
date = "04042015"
analysis.path = file.path(local.path, "analysis", date)
path.to.Rdata = file.path(local.path, "metabolome.data_032715.RData")


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
  # Creates a sample specific dataframe based on the values in sample.runs.df
  
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

save.xlsx <- function (file, ...)
{
  # Copied from http://www.r-bloggers.com/quickly-export-multiple-r-objects-to-an-excel-workbook/
  require(xlsx, quietly = TRUE)
  objects <- list(...)
  fargs <- as.list(match.call(expand.dots = TRUE))
  objnames <- as.character(fargs)[-c(1, 2)]
  nobjects <- length(objects)
  for (i in 1:nobjects) {
    if (i == 1)
      write.xlsx(objects[[i]], file, sheetName = objnames[i])
    else write.xlsx(objects[[i]], file, sheetName = objnames[i],
                    append = TRUE)
  }
  print(paste("Workbook", file, "has", nobjects, "worksheets."))
}

# Analysis


# TODO:
# Create correlation based on sample DFs
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
save.xlsx(general.fn, run.info, individual.info)

# Find the correlation for each metabolite between sample 1 and 2
correlations = diag(cor(t(as.matrix(sample.dfs[[1]])), t(as.matrix(sample.dfs[[2]]))))
hist(correlations)

# Subselect high confidence metabolites with high correlation for
# continued analysis
hf_metabolites = overlapping.metabolites[which(correlations > 0.9)]

# Create new scaled values based on sample 1 only.
all_samples = colnames(plasma[[6]])
scaled.samples = NULL
for(i in 1:length(all_samples)) {
  sample_df = NULL
  for (j in 1:length(plasma)){
    if (all_samples[i] %in% colnames(plasma[[j]])) {
      run_df = plasma[[j]][hf_metabolites, all_samples[i]]
      sample1_df = plasma[[j]][hf_metabolites, sample1]
      run_df = run_df/sample1_df
      sample_df= cbind(sample_df, run_df)
    }

  }
  rownames(sample_df) = hf_metabolites
  scaled.samples[length(scaled.samples)+1] = list(sample_df)
}
names(scaled.samples) = all_samples

cv.samples = NULL
for (i in 1:length(scaled.samples)) {
  sample_sd = apply(as.matrix(scaled.samples[[i]]), 1, sd)
  sample_mean = rowMeans(scaled.samples[[i]])
  cv_df = sample_sd/sample_mean
  cv.samples[length(cv.samples)+1] = list(cv_df)
}


