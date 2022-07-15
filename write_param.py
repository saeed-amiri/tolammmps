import itertools
import sys
import typing
import numpy as np
import pandas as pd


class Doc:
    """write LAMMPS Force Filed information
    Input:
        class object with DataFrame of Atoms, Bonds, Angles, Dihedrals
            interactions
    Output:
        LAMMPS include file
     """

class WriteParam:
    """Writting the interactions parameters for input of LAMMPS
        Input:
        atoms_df (DataFrame from PDBFILE: Pdb class)
        bonds_df, angles_df, dihedrals, impropers_df (DataFrame from
        PSFFILE Psf class)
    Output:
        A LAMMPS data file
    """
    def __init__(self, obj) -> None:
        self.obj = obj
        self.get_pairs()
        del obj
    
    def get_pairs(self) -> None:
        """find out the pairs that have somthing together"""
        # Make pair of all atom type
        type_list = self.obj.Masses_df['typ']
        pair_list = itertools.combinations_with_replacement(type_list, 2)
        
        self.write_pairs(pair_list)
    
    def write_pairs(self, pair_list: list[tuple[int, int]]) -> None:
        """write pair nonbonding interactions"""
        for pair in pair_list:
            print(pair)
