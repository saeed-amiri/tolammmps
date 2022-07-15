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
        PARAMFIEL = 'parameters.lmp'
        print(f'{self.__class__.__name__}:\n'
              f'\tWritting: "{PARAMFIEL}"\n')

        pair_list = self.mk_pairs()

        with open(PARAMFIEL, 'w') as f:
            self.write_pairs(pair_list, f)

    def mk_pairs(self) -> list[tuple[int, int]]:
        # Make pair of all atom type
        type_list = self.obj.Masses_df['typ']
        pair_list = itertools.combinations_with_replacement(type_list, 2)
        return pair_list

    def write_pairs(self,
                    pair_list: list[tuple[int, int]],
                    f: typing.TextIO) -> None:
        """write pair nonbonding interactions"""
        _df = self.obj.Masses_df.copy()
        _df.index += 1
        for pair in pair_list:
            f.write(f"# interactions between pairs\n")
            f.write(f"\n")
            f.write(f"pair_style hybrid [args...]\n")
            f.write(f"\n")
            for n, i in enumerate(pair_list):
                name_i = _df['name'][i[0]]
                name_j = _df['name'][i[1]]
                f.write(f"pair_coeff {i[0]} {i[1]} [pair_style] [args...]")
                f.write(f"#{n+1} pair: {name_i} - {name_j}\n")
            f.write(f"\n")
            del _df
