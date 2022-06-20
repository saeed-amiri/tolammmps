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


class Combine:
    """combining LAMMPS data file
    Input:
        list of the name of input files
    """

    def __init__(self, file_list: list[str]) -> None:
        self.f_list = file_list

    def mk_lmp_df(self) -> None:
        """making DataFrames from all the inputs"""
        self.set_df_lists()
        self.get_data()
        self.update_atoms_df()

    def set_df_lists(self) -> None:
        """Set lists to append DataFrame in it"""
        self.l_atoms: list[pd.DataFrame] = []
        self.l_bonds: list[pd.DataFrame] = []
        self.l_angles: list[pd.DataFrame] = []
        self.l_dihedrals: list[pd.DataFrame] = []
        self.l_headers: dict[str, mlmp.Header] = dict()

    def get_data(self) -> None:
        """loop over all the files and make several DataFrame"""
        for f in self.f_list:
            print(f)
            o = mlmp.Header(f)
            self.l_headers[f] = o
            atoms = mlmp.Body(o.Names, f)
            atoms.read_body()
            self.l_atoms.append(atoms.Atoms_df)
            self.l_bonds.append(atoms.Bonds_df)
            self.l_angles.append(atoms.Angles_df)
            self.l_dihedrals.append(atoms.Dihedrals_df)

    def append_atoms(self) -> None:
        """append atoms DataFrame with updataed id and type"""
        for df in self.l_atoms:
            print(df)

    def update_atoms_df(self) -> None:
        """Update the type and index of the atoms with number of next
        atoms in the next DataFrame
        """
        



INFILE = sys.argv[1:]
system = Combine(INFILE)
system.mk_lmp_df()
print(next(iter(system.l_headers)))
