import pandas as pd
import scipy.stats as stats
import helpers
from change_detector import ZScoreDetector, OnlineSimulator


def calculate_zscores(series, sample=None,
                      include_sample=False):
    '''
    Calculate z scores for an entire series
    or for a particular value. The sample can
    be part of population or not.

    :param series: numpy.array_like values
    :param sample: returns zscores for only that sample.
    :param include_sample: include sample in calculation of mean and std
    :returns: zscores
    '''
    if sample:
        value = series.ix[sample]
        if not include_sample:
            series = series.ix[set(series.index) - set([sample])]
        return helpers.zscore(value, series.mean(), series.std(ddof=1))
    else:
        return stats.zscore(series)


def outliers(df, z_threshold=2,
             sample=None, include_sample=False):
    '''Returns the classifications for each sample,
     for which it is above a z-score threshold.

    :param z_threshold: absolute value z threshold
    :param sample: returns zscores for only that sample.
    :param include_sample: include sample in calculation of mean and std
    :returns
    Returns classifications
    '''

    zscores = df.apply(calculate_zscores,
                       args=(sample, include_sample))
    if sample:
        zscores = pd.DataFrame(zscores).T
        zscores.index = [sample]
    zscores = zscores.fillna(0).applymap(abs)
    outliers = zscores > z_threshold
    # Hack, because multindex columns apply in pandas is broken.
    # outliers = outliers.apply(lambda x: ([y for y in x.index if x[y]]), axis=0)
    outlier_columns = outliers.columns
    outliers.columns = range(outliers.shape[1])
    outliers = outliers.apply(lambda x: ([y for y in x.index if x[y]]), axis=0)
    outliers.index = outlier_columns
    outliers = pd.DataFrame(outliers, columns=['outliers'])
    return outliers


def regress(df,):
    '''
    Returns linear regression p value for each phenotype column in the df
    '''
    days = df.index.map(lambda x: x - df.index.min())
    days = days.days
    # pvalues = df.apply(lambda x: stats.linregress(days, x).pvalue)
    # 4th entry is pvalue
    pvalues = pd.DataFrame(df.apply(lambda x: stats.linregress(days, x)[3]),
                           columns=['linear-pvalue'])
    return pvalues


def changepoints(df, window_size=5, threshold=1.67):
    '''
    Returns change points for along every classification
    '''
    def detect(series):
        detector = ZScoreDetector(window_size, threshold)
        simulator = OnlineSimulator(detector, series)
        change = simulator.run(plot=False)
        if change:
            changepoint = simulator.residuals_history.itervalues().next()
            changepoint = len(changepoint) - window_size
            return series.index[changepoint]
        return None

    changepoints = pd.DataFrame(df.apply(detect), columns=['changepoints'])
    return changepoints
