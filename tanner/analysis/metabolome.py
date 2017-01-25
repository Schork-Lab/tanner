"""
Metabolomics analysis functions
"""
import pandas as pd
import tanner.stats.bayesian as bay
import pymc3 as pm
import os
import sys

def load_df(fn, drop_missing=True):
    df = pd.read_table(fn, index_col=[0, 1], parse_dates=True)
    if drop_missing:
        df = df.dropna(axis=1)
    return df

def parse_metabolite(df, metabolite):
    """
    Parse a metabolite dataframe into constituent parts.
    """
    df = df[[metabolite]]
    df.columns = ['metabolite']
    min_day = df.index.get_level_values(0).min()
    df.loc[:,'day'] = df.index.get_level_values(0).map(lambda x: (x-min_day).days)
    df.loc[:,'run'] = df.index.get_level_values(1)
    median_metabolites = df.groupby('run').median()['metabolite']
    df.loc[:, 'median_scaled'] = df.apply(lambda x: x['metabolite']/median_metabolites[x['run']], axis=1)
    df = df.sort_values(by=['run', 'day'])
    return df

def bayesian_fit(fn, metabolite_col, out_dir, min_run=4, n_chain=50000, pooled=False):

    df = load_df(fn, drop_missing=False)
    metabolite = df.columns[metabolite_col]
    print('Running Bayesian Full Model for {}'.format(metabolite))
    metabolite_df = parse_metabolite(df, metabolite).dropna()
    model = bay.Linear(variational=False, n_chain=n_chain, pooled=pooled)
    parsed_dict = bay.parse_df(metabolite_df, min_run)
    model.run(**parsed_dict)
    summary = model.summary()
    pooled_phrase = "pooled" if pooled else "not_pooled"
    out_fn = os.path.join(out_dir, '{}.{}.info.tsv'.format(metabolite, pooled_phrase))
    pd.DataFrame.from_dict(parsed_dict, orient='index').to_csv(out_fn, sep='\t')
    summary.to_csv(out_fn.replace('info','bay_summary'), sep='\t')
    return model.trace

def simple_fit(fn, metabolite_col, out_dir):
    df = load_df(fn, drop_missing=False)
    metabolite = df.columns[metabolite_col]
    metabolite_df = parse_metabolite(df, metabolite).dropna()

    for level in (8, 9):
        parsed_dict = {'time_values': metabolite_df.xs(level, level=1)['day'].values,
                       'measured_levels': metabolite_df.xs(level, level=1)['metabolite'].values}

        model = bay.SimpleLinear(variational=False)
        model.run(**parsed_dict)
        linear_summary = model.summary()
        
        model = bay.SimpleLinearNoScale(variational=False)
        model.run(**parsed_dict)
        noscale_linear_summary = model.summary()

        out_fn = os.path.join(out_dir, '{}.{}.simple.info.tsv'.format(metabolite, level))
        pd.DataFrame.from_dict(parsed_dict, orient='index').to_csv(out_fn, sep='\t')
        linear_summary.to_csv(out_fn.replace('info','linear_summary'), sep='\t')
        noscale_linear_summary.to_csv(out_fn.replace('info','noscale_linear_summary'), sep='\t')

    return model.trace

def __main__():
    if len(sys.argv) < 4:
        print("Usage: metabolome-fit not-pooled fn metabolite_col out_dir")
        print("Usage: metabolite-fit pooled fn metabolite_col out_dir")
        print("Usage: metabolite-fit SIMPLE fn metabolite_col out_dir")
    else:
        if sys.argv[1] == 'not-pooled':
            print("Running Bayesian Non-pooled for {}".format("\t".join(sys.argv[1:])))
            bayesian_fit(sys.argv[2], int(sys.argv[3]), sys.argv[4])
        elif sys.argv[1] == 'pooled':
            print("Running Bayesian Pooled for {}".format("\t".join(sys.argv[1:])))
            bayesian_fit(sys.argv[2], int(sys.argv[3]), sys.argv[4], pooled=True)
        else:    
            print("Running Simple Bayesian fit for {}".format("\t".join(sys.argv[1:])))
            simple_fit(sys.argv[2], int(sys.argv[3]), sys.argv[4])