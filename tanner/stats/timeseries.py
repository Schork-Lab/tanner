import pandas as pd
import numpy as np
import scipy.stats as stats
import statsmodels.api as sm
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
import tanner.stats.helpers as helpers

class Statistics():
    def __init__(self, ts, time=None, window_size=3, threshold=1.67,
                 min_time=5, name=None):
        self.ts = ts
        if time is None:
            time = list(range(len(ts)))
        self.time = time
        self.window_size = window_size
        self.threshold = threshold
        self.min_time = min_time
        self.name = name
        self.cp_, self.outliers_, self.regress_ = [], [], [] 
        return

    def df(self):
        if self.cp_ == []:
            self._calc_statistics()
        df = pd.DataFrame.from_dict({'outliers': self.outliers_,
                                     'changepoint': self.cp_,
                                     'regression': self.regress_},
                                     orient='index')
        df.columns = self.time
        if self.name is not None:
            df.index = pd.MultiIndex.from_tuples([(self.name, index) for index in df.index])

        return df

    def run(self):
        self._calc_statistics()

    def _calc_statistics(self):
        """
        override this function
        """


class TStatistics(Statistics):
    def __init__(self, *args, **kwargs):
        super(TStatistics, self).__init__(*args, **kwargs)

    def _calc_statistics(self):
        self.cp_, self.outliers_, self.regress_ = [], [], [] 
        for i, time in enumerate(self.time):
            if i < self.min_time:
                self.cp_.append(np.nan)
                self.outliers_.append(np.nan)
                self.regress_.append(np.nan)
                continue 

            ts = self.ts[:i+1]

            # mean difference changepoint
            # TODO: g_std = g_std/sqrt(n)? sample std instead of population
            window = ts[-self.window_size:]
            g_mean, g_std = np.mean(ts), np.std(ts)
            w_mean, w_std = np.mean(window), np.std(window)    
            SE = g_std / np.sqrt(self.window_size)
            cp_tstat = (w_mean - g_mean) / SE
            cp_pvalue = 2*stats.t.sf(np.abs(cp_tstat), df=i) #i-self.window_size)?
            self.cp_.append((i-self.window_size, cp_pvalue))

            # outlier
            mean, sd = np.mean(ts[:-1]), np.std(ts[:-1])
            outlier_tstat = (ts[-1] - mean) / (sd/np.sqrt(i-1))
            outlier_pvalue = 2*stats.t.sf(np.abs(outlier_tstat), df=i-1)
            self.outliers_.append(outlier_pvalue)

            # regression
            days = self.time[:i+1]
            slope, intercept, r_value, p_value, std_err = stats.linregress(days, ts)
            self.regress_.append(p_value)


class RobustStatistics(Statistics):
    def __init__(self, RLM=True, *args, **kwargs):
        if RLM:
            self.linear_estimator = sm.RLM
        else:
            self.linear_estimator = sm.OLS
        super(RobustStatistics, self).__init__(*args, **kwargs)
    
    def _calc_statistics(self):
        for i, time in enumerate(self.time):
            if i < self.min_time:
                self.cp_.append(np.nan)
                self.outliers_.append(np.nan)
                self.regress_.append(np.nan)
                continue 

            ts = self.ts[:i+1]
            days = self.time[:i+1]

            # Changepoint
            ecp, pval = self._estimate_cp_pval(ts)
            if ecp == i+1:
                self.cp_.append(np.nan)
            else:
                self.cp_.append((ecp, pval))

            # outlier
            # http://stackoverflow.com/questions/22354094/pythonic-way-of-detecting-outliers-in-one-dimensional-observation-data
            # modified_z_score = self._mad_outlier(ts)
            # self.outliers_.append(2*stats.norm.sf(np.abs(modified_z_score[-1])))
            # studentized outlier
            pvalue = self._studentized_outlier(ts, days)
            self.outliers_.append(pvalue)

            # regression
            
            pvalue = self._robust_regress(ts, days)
            self.regress_.append(pvalue)

    def _robust_regress(self, ts, days):
        design = sm.add_constant(days)
        fit = self.linear_estimator(ts, design).fit()
        return fit.pvalues[-1]

    def _studentized_outlier(self, ts, days):
        design = sm.add_constant(days)
        fit = sm.OLS(ts, design).fit()
        outliers = fit.outlier_test()
        return outliers.ix[-1]['unadj_p']

    def _mad_outlier(self, ts):
        median = np.median(ts)
        diff = np.abs(ts - median)
        med_abs_deviation = np.median(diff)
        modified_z_score = 0.6745 * diff / med_abs_deviation
        return modified_z_score

    def _estimate_cp_pval(self, ts, method="mean"):
        """ https://github.com/viveksck/changepoint/blob/master/changepoint/rchangepoint.py """
        """ Estimate changepoints in a time series by using R. """

        """ 
            ts: time series
            method: look for a single changepoint in 'mean' , 'var', 'mean and var'
            Returns: returns index of the changepoint, and pvalue. Here pvalue = 1
                     means statistically significant
        """
        robjects.r("library(changepoint)")
        method_map = {
            "mean": "cpt.mean({}, class=FALSE)",
            "var": "cpt.var({}, class=FALSE)",
            "meanvar": "cpt.meanvar({}, class=FALSE)",
        }
        mt = robjects.FloatVector(ts)
        robjects.globalenv["mt"] = mt
        cmd = method_map[method].format("mt")
        robjects.globalenv["mycpt"] = robjects.r(cmd)
        ecp_pval = robjects.r("mycpt")
        ecp = ecp_pval[0]
        pval = ecp_pval[1]
        return ecp, pval


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


# def changepoints(df, window_size=5, threshold=1.67):
#     '''
#     Returns change points for along every classification
#     '''
#     def detect(series):
#         detector = ZScoreDetector(window_size, threshold)
#         simulator = OnlineSimulator(detector, series)
#         change = simulator.run(plot=False)
#         if change:
#             changepoint = list(simulator.residuals_history.values())[0]
#             changepoint = len(changepoint) - window_size
#             return series.index[changepoint]
#         return None

#     changepoints = pd.DataFrame([df.apply(detect)]).dropna(axis=1).T
#     changepoints.columns = ['changepoints']
#     #changepoints = pd.DataFrame(df.apply(detect), columns=['changepoints'])
#     return changepoints
