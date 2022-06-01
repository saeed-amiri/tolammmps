import pandas as pd

"""
Reading data file from Ole Nickel, received on Jun 01, 2022
Data contains many atoms,
The atoms with types 1 to 4 are the ones that belong to the SiO2
input:
    - LAMMPS data file
output:
    - LAMMPS data file

I'm going to use fast writing, with setting stderr to file
"""

class ReadData:
    """reading file
    Ignore the header and then read the ATOM cards
    """
    pass


INFILE = 'merged.data'

