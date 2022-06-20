import os
import sys
import pandas as pd
import numpy as np

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