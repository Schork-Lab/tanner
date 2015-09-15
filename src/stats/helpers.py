from numpy import log10

zscore = lambda value, mean, std: (mean - value)/std

def qqvalues(pvalues):
    x = [-log10(float(i)/len(pvalues)) 
         for i in range(1, len(pvalues)+1)]
    x.sort()
    y = pvalues.map(lambda x: -log10(x))
    y.sort()
    return x, y