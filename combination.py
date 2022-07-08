import sys
import os
import re
import pandas as pd

class Doc:
    """This code combines data files and returns a superstructure in
    LAMMPS full atom style.   The input file is <system>.struct file
    which contains the names and path of the  files and also a matrix
    that shows the order of the superstructure:

   ex. of an input file:
    
    First, it should have a symbol for each file and its path.
    The symbol is better to be an upper letter (for now, it does not
    recognize lower/upper case of the letter):

        # D=decane.data
        # W=water.data
        # S=sio2.data
    Then it should have a matrix that shows how you want to build your
    superstructure:
        DWD
        _S_
        DWD
    which stack the blocks as: "decane water decane" most top level,
    SiO2 will stack in the second layer and "decane water decane" in
    the lowest layer.

    Few reserved char:
    | (pipe char): which means the upper block continues vertically to
    here. 
    _ (underline char): means that the previous structure continues
    here.
    - (dash): empty space

    Jul 08 2022
    Saeed
    """

class Structure:
    """read the struct file"""
    def __init__(self) -> None:
        pass