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
        self.__header: str  # str to print on the parameter file sections
        self.__header = ''.join(['#']*79)
        pair_list = self.mk_pairs()
        self.mk_bond_pair()

        with open(PARAMFIEL, 'w') as f:
            self.write_pairs(pair_list, f)

    def mk_pairs(self) -> list[tuple[int, int]]:
        # Make pair of all atom type
        type_list = self.obj.Masses_df['typ']
        pair_list = itertools.combinations_with_replacement(type_list, 2)
        return pair_list

    def mk_bond_pair(self) -> pd.DataFrame:
        """make a DataFrame for type and name of bond interactions"""
        # Clairfy bond interactions
        _df: pd.DataFrame  # Temperary Dataframe to group Bonds_df
        bond_list: list[int]  # List of bond types
        ai_type: list[int]  # List of atom types
        aj_type: list[int]  # List of atom types
        ai_name: list[str]  # List of atoms name
        aj_name: list[str]  # List of atoms name
        i_id: int  # To temp save the atoms type
        j_id: int  # To temp save the atoms type
        _df_bpair: pd.DataFrame  # To return bond pair information

        _df = self.obj.Bonds_df.copy()
        _df = _df.groupby(by='typ').min()
        _df = _df.reset_index()  # To make "typ" as a column
        bond_list = _df['typ']
        _df_bpair = pd.DataFrame(bond_list)
        _df_bpair = _df_bpair.rename(columns={'typ': 'bond_typ'})
        ai_type = []
        aj_type = []
        ai_name = []
        aj_name = []
        for ai, aj in zip(_df['ai'], _df['aj']):
            i_id = self.obj.Atoms_df.iloc[ai-1]['typ']
            j_id = self.obj.Atoms_df.iloc[aj-1]['typ']
            ai_type.append(i_id)
            aj_type.append(j_id)
            ai_name.append(self.obj.Masses_df.iloc[i_id-1]['name'])
            aj_name.append(self.obj.Masses_df.iloc[j_id-1]['name'])
        _df_bpair['ai_typ'] = ai_type
        _df_bpair['aj_typ'] = aj_type
        _df_bpair['ai_name'] = ai_name
        _df_bpair['aj_name'] = aj_name
        del _df, ai_type, aj_type,  ai_name,  aj_name
        return _df_bpair

    def write_pairs(self,
                    pair_list: list[tuple[int, int]],
                    f: typing.TextIO) -> None:
        """write pair nonbonding interactions"""
        _df = self.obj.Masses_df.copy()
        _df.index += 1
        for pair in pair_list:
            f.write(f'{self.__header}\n')
            f.write(f'# interactions between pairs\n')
            f.write(f'\n')
            f.write(f'pair_style hybrid [args...]\n')
            f.write(f'\n')
            for i, pair in enumerate(pair_list):
                name_i = _df['name'][pair[0]]
                name_j = _df['name'][pair[1]]
                f.write(f'pair_coeff {pair[0]} {pair[1]}'
                        f' [pair_style] [args...]')
                f.write(f' # {i+1} pair: {name_i} - {name_j}\n')
            f.write(f"\n")
        del _df
