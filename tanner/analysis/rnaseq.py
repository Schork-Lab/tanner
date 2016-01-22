'''
Analysis code related to the newest RNASeq samples
'''
import os
import pandas as pd
import helpers


def load_deseq(directory, deseq_type="normalized", individual=None):
    """
    Load DESeq log normalized data for further analysis. The choice
    of log transformed data was made based on reading the following discussions

    http://seqanswers.com/forums/showthread.php?t=41421
    https://support.bioconductor.org/p/61884/#61895

    In general, it's to stabilize the variance, so that we do not
    identify outliers in downstream time series analysis that are
    only there due to biases associated with lower-abundance transcripts

    """
    annotated_count_fn = os.path.join(directory, "annotated_combined.counts")
    count_df = pd.read_table(annotated_count_fn, sep="\t")
    if deseq_type == "normalized":
        deseq_file = os.path.join(directory, "deseq.counts")
    elif deseq_type == "rlog":
        deseq_file = os.path.join(directory, "deseq.regularized.log.counts")
    elif deseq_type == "vstab":
        deseq_file = os.path.join(directory, "deseq.variance.stablized.counts")
    else:
        raise ValueError('Unrecognized deseq file type')
    deseq_df = pd.read_table(deseq_file, sep="\t")
    deseq_df['symbol'] = count_df['symbol'].values
    deseq_df = deseq_df.drop_duplicates(subset='symbol', take_last=True)
    deseq_df = deseq_df.set_index('symbol')
    deseq_df = deseq_df.T
    samples = deseq_df.index.map(lambda x: x.split('_')[0])
    dates = pd.to_datetime(deseq_df.index.map(lambda x: x.split('_')[-1]))
    deseq_df.index = pd.MultiIndex.from_tuples(zip(samples, dates))
    if individual:
        return deseq_df.ix[individual, :]
    else:
        return deseq_df


def process_timeseries(df, analysis_path, pvalue=0.05):
    '''
    Process the basic time series analysis including
    outlier detection, linear trends, and changepoints.
    '''
    helpers.process_timeseries(df, analysis_path, pvalue)
