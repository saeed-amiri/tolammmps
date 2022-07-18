import re
import os
import sys
import pandas as pd


class Doc:
    """read the 'param.data': parameter file to write LAMMPS
    parameter.lmp file
        Input:
            param.data, which should contain:
                pair style
                style coefficients
                bond, angle, dihedral styles, and parameters
        Out:
            pd.Dataframe contains all the parameters of interactions.
    """


class ReadParameter:
    """read parmater file"""
