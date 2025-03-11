import numpy as np

def get_cdf(data, type = 'cdf'):
    vals = np.unique(data)
    cnts = np.array([sum(item == data) for item in vals])
    csum = np.cumsum(cnts)
    cdf  = csum/csum[-1]
    pdf  = cnts/sum(cnts)

    df = cdf if type == 'cdf' else pdf

    return vals, df
    
