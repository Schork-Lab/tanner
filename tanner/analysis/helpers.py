import os
import scipy.stats as stats
import tanner.stats.timeseries as ts
import tanner.visual.timeseries as vts
import pandas as pd

def create_timeseries_figures(data_df, output_folder, pvalue=0.05):
    '''
    Process the basic time series analysis including
    outlier detection, linear trends, and changepoints.

    Save all interesting plots in the output_folder
    '''
    # Zscore cutoffs:
    # pvalue = pvalue/len(data_df.columns)
    zscore = abs(stats.norm.ppf(pvalue))

    # Outliers
    outliers = ts.outliers(data_df, z_threshold=zscore)
    for sample, classification in outliers.iteritems():
        if len(classification) > 0:
            date = str(sample).split(' ')[0]
            file_prefix = os.path.join(output_folder, date+'_outliers')
            vts.default(data_df[classification], fn_prefix=file_prefix,
                        main_title=sample, marksample=sample)
    outliers_df = pd.DataFrame.from_dict(outliers)


    # Linear Regression:
    pvalues = ts.regress(data_df)
    file_prefix = os.path.join(output_folder, 'lineartrends')
    significant = pvalues[pvalues < pvalue].index
    if len(significant) > 0:
        vts.default(data_df[significant],
                    fn_prefix=file_prefix,
                    bestfit=True)

    # Changepoints:
    file_prefix = os.path.join(output_folder, 'changepoints')
    changepoints = ts.changepoints(data_df, threshold=zscore)
    changepoints = changepoints[changepoints > 0]
    if len(changepoints) > 0:
        vts.default(data_df[changepoints.index],
                    fn_prefix=file_prefix,
                    main_title='Changepoint',
                    changepoints=changepoints.values)


def create_timeseries_df(data_df, pvalue):
    '''
    '''
    def time_to_str(x):
        try:
            return x.strftime('%Y-%m-%d')
        except ValueError:
            return 'None'
    zscore = abs(stats.norm.ppf(pvalue))
    zscores = data_df.apply(ts.calculate_zscores)
    outliers = ts.outliers(data_df, zscore)
    pvalues = ts.regress(data_df)
    trends = pvalues < pvalue
    trends.columns = ['linear-trend']
    if zscore > 2:
        zscore = 2
    changepoints = ts.changepoints(data_df, threshold=zscore)
    zscores = zscores.fillna(0).apply(lambda x: ",".join([str(y) for y in x]))
    zscores = pd.DataFrame(zscores, columns=['zscores'])
    outliers = outliers.applymap(lambda x: ",".join([time_to_str(y) for y in x]))
    outliers[outliers == ''] = 'None'
    changepoints = changepoints.applymap(time_to_str).fillna('None')
    ts_df = pd.concat([zscores, outliers, changepoints, trends, pvalues], axis=1)
    return ts_df
