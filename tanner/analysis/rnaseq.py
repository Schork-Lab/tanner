'''
Analysis code related to the newest RNASeq samples
'''
import os
import pandas as pd
import helpers


def load_deseq(directory, individual=None):
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
    rlog_fn = os.path.join(directory, "deseq.regularized.log.counts")
    rlog_df = pd.read_table(rlog_fn, sep="\t")
    rlog_df['symbol'] = count_df['symbol'].values
    rlog_df = rlog_df.drop_duplicates(subset='symbol', take_last=True)
    rlog_df = rlog_df.set_index('symbol')
    rlog_df = rlog_df.T
    rlog_samples = rlog_df.index.map(lambda x: x.split('_')[0])
    rlog_dates = pd.to_datetime(rlog_df.index.map(lambda x: x.split('_')[-1]))
    rlog_df.index = pd.MultiIndex.from_tuples(zip(rlog_samples, rlog_dates))
    if individual:
        return rlog_df.ix[individual, :]
    else:
        return rlog_df


def process_timeseries(df, analysis_path, pvalue=0.05):
    '''
    Process the basic time series analysis including
    outlier detection, linear trends, and changepoints.
    '''
    helpers.process_timeseries(df, analysis_path, pvalue)
