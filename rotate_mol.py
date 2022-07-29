import os
import sys
import pandas as pd
import numpy as np
import read_lmp_data as rlmp
import write_lmp as wlmp


class Doc:
    """rotate molecules
    Rotate molecules around thier axis randomly
    Input:
        LAMMPS data file
    Output:
        LAMMPS data file
    """


class Rotate:
    """get the data file and rotate along axis"""
    def __init__(self) -> None:
        pass


if __name__ == '__main__':
    fname = sys.argv[1]
    data = rlmp.ReadData(fname)
    print(data.Atoms_df)
