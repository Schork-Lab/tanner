import os
import scipy.stats as stats
import tanner.stats.timeseries as ts
import tanner.visual.timeseries as vts


def process_timeseries(data_df, output_folder, pvalue_threshold=0.05):
    '''
    Process the basic time series analysis including
    outlier detection, linear trends, and changepoints.

    Save all interesting plots in the output_folder
    '''
    # Zscore cutoffs:
    pvalue = pvalue_threshold/len(data_df.columns)
    zscore = abs(stats.norm.ppf(pvalue))

    # Outliers
    outliers = ts.outliers(data_df, z_threshold=zscore)
    for sample, classification in outliers.iteritems():
        if len(classification) > 0:
            date = str(sample).split(' ')[0]
            file_prefix = os.path.join(output_folder, date+'_outliers')
            vts.default(data_df[classification], fn_prefix=file_prefix,
                        main_title=sample, marksample=sample)

    print pvalue

    # Linear Regression:
    pvalues = ts.regress(data_df)
    print pvalues.head().values
    file_prefix = os.path.join(output_folder, 'lineartrend')
    significant = pvalues[pvalues < pvalue].index

    print sum(pvalues < pvalue)
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
