import os, sys, re
import pandas as pd

class DOC:
    """"
    Reading the AMBER data file for SiO2 slab and converting to LAMMPS data file
    """

class TOP:
    """
    reading top file
    it has 40 FLAG (cards) with different formats in FORTAN style written in the next line
    each flag:
    %FLAG POINTERS                                                                  
    %FORMAT(10I8)                                                                   
    6702       7    5350       0       0       0       0       0       0       0
    8468    3062       0       0       0       3       0       0      10       1
       0       0       0       0       0       0       0       1       8       0
       0

    The cards are space-parsed in some cases, i.e. names of atoms, mols, residues
    """
    def __init__(self) -> None:
        pass

    def read_file(self):
        """
        read line by line and subtract the data from them
        """
        with open 


if __name__=="__main__":
    TOPFILE = "test.top"