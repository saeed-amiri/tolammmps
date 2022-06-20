import os
import sys
import typing
import pandas as pd
import numpy as np
import read_lmp_data as mlmp  # My lammps

class Doc:
    """combining LAMMPS data file to prepare interface
    Input:
        LAMMPS data file[s]
    Output:
        LAMMPS data file
    
    The order of data is based on the input order.
    For now, all the data only stack on top of each other in the
    z-direction.

    Usage:
        combine_data.py data1 data2 data3 ...
    
    """


INFILE = 'system.data'
OUTFILE = 'silica_ole_slic.data'
o = mlmp.Header(INFILE)
atoms = mlmp.Body(o.Names, INFILE)
print(o.Names)
atoms.read_body()
