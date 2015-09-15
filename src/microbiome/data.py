import pandas as pd


def load_aggregated(fn, num_of_samples=11):
    '''
    Loads aggregated files
    Sets indices as datetime.
    '''
    assert(fn.endswith(".txt"))
    assert(num_of_samples > 0)
    df = pd.read_table(fn,)
    df = df[[column for column in df.columns
            if (column.find('ID') == -1) and
            (column.find('No_families') == -1)]]
    df.set_index(list(df.columns[:-num_of_samples]), drop=True, inplace=True)
    df = df.T
    df.index = pd.to_datetime(df.index.map(lambda x: x.split('_')[-1]))
    return df
