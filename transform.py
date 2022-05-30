import pandas as pd
import itertools
import numpy as np

"""
A script for transform parameters table to lammps pair interaction
#v Si OB OH H OM
t 1 2 3 4 5
q 1.600 -0.800 -0.800 0.400 -1.000
s 1.600 1.762 1.650 1.000 1.650
e 0.300 0.261 0.150 0.021 0.150
"""

def read_param(filename) -> list:
    params= []
    header = []
    with open (filename, 'r') as f:
        while True:
            line = f.readline()
            if line.startswith('#'): header.append(line.strip())
            else: params.append(line.strip())
            if not line: break
    return header, params

def mk_df(header, params) -> pd.DataFrame:
    df = pd.DataFrame(params, columns=header)
    for i, c in enumerate(header):
        if i!=0: df = df.astype({c:float})
    df  = df.set_index("#v")
    return df

def mk_table(df) -> pd.DataFrame:
    for i, j in itertools.combinations_with_replacement(df.keys(), 2):
        # print(i, j)
        epsilon, sigma = mk_pairs(df[i], df[j])
        print(i, j, epsilon, sigma)

def mk_pairs(ai, aj) -> float:
    epsilon = np.sqrt(ai.e * aj.e)
    sigma = (ai.s + aj.s )*0.5
    return epsilon, sigma

# get header
header, params = read_param('anke.param')
header = header[0].split(' ')
# get parameters
params = [item.split() for item in params]
params = [item for item in params if item]
# make dataframe
df = mk_df(header, params)
# make lammps pair interaction table
mk_table(df)