#Author: Kunal Bhutani
#Date: March 31s, 2015
#This script does basic QC for the metabolome data. 


path.to.Rdata = "~/tscc_projects/tanner/data/family3/metabolome/metabolome.data_032715.RData"
load(path.to.Rdata)
GetAllSamples <- function() {
  samples = NULL
  for (i in 1:length(plasma)) {
    run.samples = names(plasma[[i]])[grepl("X001", names(plasma[[i]]))]
    print(run.samples)
    samples = union(samples, run.samples)
  }
  samples
}

samples = GetAllSamples()

sample1 = 'X001_002_Metabolites_Plasma_JCVI.00001_04.16.14'
sample2 = 'X001_002_Metabolites_Plasma_JCVI.00002_04.16.14'
samples = c(sample1, sample2)


FindSampleRuns <- function(samples) {
  # For each sample, create a list of runs that the sample is a part of
  relevant.runs = NULL
  for (i in 1:length(plasma)) {
    m = match(samples, colnames(plasma[[i]]))
    if (sum(!is.na(m)) == length(samples)){
      relevant.runs = c(relevant.runs, i)
    }
   }
  relevant.runs
}

FindOverlappingMetabolites <- function(runs) {
  # Finds a set of metabolites that are found in each of the runs
  overlapping.metabolites = rownames(plasma[[runs[1]]])
  for (i in 2:length(runs)){
    overlapping.metabolites = intersect(overlapping.metabolites, 
                                        rownames(plasma[[runs[i]]]))
  }
  overlapping.metabolites
}


CreateSampleDfs <- function(samples, runs, overlapping.metabolites) {
  # Create sample specific dataframes
  sample.dfs = NULL
  for(i in 1:length(samples)) {
    sample.df = data.frame(plasma[[runs[1]]][overlapping.metabolites, samples[i]])
    colnames(sample.df) = c(paste("Run", runs[1]))
    for (j in 2:length(runs)){
      run.df = data.frame(plasma[[runs[j]]][overlapping.metabolites, samples[i]])
      colnames(run.df) = c(paste("Run", runs[j]))
      sample.df= cbind(sample.df, run.df)
    }
    rownames(sample.df) = overlapping.metabolites
    sample.dfs[length(sample.dfs)+1] = list(sample.df)
  }
  names(sample.dfs) = samples
  sample.dfs
}

relevant.runs = FindSampleRuns(samples)
overlapping.metabolites = FindOverlappingMetabolites(relevant.runs)
sample.dfs = CreateSampleDfs(samples, relevant.runs, overlapping.metabolites)


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


