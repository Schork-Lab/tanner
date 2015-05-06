#Author: Kunal Bhutani
#Date: March 31st, 2015
#Last Updated: May 5th, 2015
# Does the initial loading of all metabolomics data.

# Libraries
library(xlsx)
library(yaml)

# Configuring Paths
config = yaml.load_file('../config.yaml')
file.dir = dirname(parent.frame(2)$ofile)
setwd(file.dir)
path.to.Rdata = "data/family3/metabolome/metabolome.data_032715.RData"
path.to.Rdata = file.path(config$paths$local, path.to.Rdata)


# Loading in Data
load(path.to.Rdata)
rownames(metab.info) = metab.info[,'CHEMICAL.ID']

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
    sample.runs = which(!is.na(sample.runs.df[sample, ]))
    sample.df = data.frame(lapply(sample.runs, function(x) plasma[[x]][metabolites,
                                                                          sample]))
    rownames(sample.df) = metabolites
    colnames(sample.df) = sample.runs
    sample.df
  }
  sample.dfs = lapply(rownames(sample.runs.df), ProcessSample)
  names(sample.dfs) = rownames(sample.runs.df)
  sample.dfs
}

CreateMetaboliteDfs <- function(sample.runs.df, sample.dfs,
                                metabolites) {
  # Rows: Each different sample
  # Sample order is determined by going across the sample.runs.df one row at a time
  # Columns: Metabolites
  # Special Columns: Identities/Factors for Run #
  # Special Columns: Identities/Factors for Sample
  
  # Get number of samples and number of runs from sample.runs.df
  num.samples = dim(sample.runs.df)[0]
  num.runs =dim(sample.runs.df)[1]
  num.total = sum(!is.na(as.vector(t(sample.runs.df))))
  
  # Since the metabolite df is melted down to 2-d, the ordering is determined by looking
  # across runs for each sample. Then down each row.
  sample.ids.order = unlist(lapply(rownames(sample.runs.df), 
                                   function (x) { rep(x, length(which(!is.na(sample.runs.df[x,]))))}
                                   ))
  
  run.ids.order = unlist(lapply(as.data.frame(t(sample.runs.df)), function (x) { x[which(!is.na(x))] }))
  run.ids.order = colnames(sample.runs.df)[sample.order]
  
  # Create column factors for samples and columns for mixed modeling and other applications.
  sample.cols = as.factor(sample.ids.order)
  run.cols = as.factor(run.ids.order)

  
# The commented out is a different method of creating above but it also creates a different
# sample and run column for each sample and run instead of a single factor. This might be useful in other
# aplications.
#   run.cols = sapply(colnames(sample.runs.df),
#                     function (x) {
#                       run = rep(0, num.total)
#                       run[which(run.ids.order == x)] = 1
#                       as.factor(run)
#                     })
#   colnames(run.cols) = sapply(colnames(sample.runs.df), function (x) {paste("Run", x, sep="X")})
#   
#   sample.cols = sapply(rownames(sample.runs.df),
#                          function (x) {
#                            sample = rep(0, num.total)
#                            sample[grepl(x, sample.ids.order)] = 1
#                            as.factor(sample)
#                          })

  # Create a column for each metabolite based on the previously determined order for rows
  ProcessMetabolite = function(metabolite) {
    values = sapply(1:num.total,
                    function(i) { sample.dfs[[sample.ids.order[i]]][metabolite, run.ids.order[i]]})
    values
  }
  metabolite.cols = as.data.frame(sapply(metabolites, ProcessMetabolite))
  colnames(metabolite.cols) = metabolites

  # Combine all the information into a single 
  combined.df = as.data.frame(cbind(run.cols, sample.cols, metabolite.cols))
  colnames(combined.df) = unlist(lapply(colnames(metabolite.df), function(x) {paste("X", x, sep="")}))
  combined.df
  
}

# Loading default parameters
samples = GetAllSamples()
sample.runs.df = CreateSampleRunsDf(samples)
overlapping.metabolites = FindOverlappingMetabolites()
sample.dfs = CreateSampleDfs(sample.runs.df, overlapping.metabolites)
metabolite.df = CreateMetaboliteDfs(sample.runs.df, sample.dfs, overlapping.metabolites)

