"""
Metabolomics analysis functions
"""
import pandas as pd
import tanner.stats.bayesian as bay
import pymc3 as pm
import os

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

def bayesian_fit(fn, metabolite_col, out_dir, min_run=4):

    df = load_df(fn)
    metabolite = df.columns[metabolite_col]
    metabolite_df = parse_metabolite(df, metabolite)
    model = bay.Linear(variational=False)
    model.run(**bay.parse_df(metabolite_df, min_run))
    summary = model.summary()
    out_fn = os.path.join(out_dir, '{}.bay_summary.tsv'.format(metabolite))
    summary.to_csv(out_fn, sep='\t')

    return model.trace