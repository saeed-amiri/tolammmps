import re
import typing
import itertools
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
    def __init__(self, obj, bs) -> None:
        self.obj = obj
        self.param = bs.param
        self.get_pairs()
        del obj, bs

    def mk_parameters(self) -> None:
        """make dataframe from updated parameter file"""
        df: pd.DataFrame  # Temporary dataframe
        param_df: pd.DataFrame  # The main dataframe to return
        df_list: list[pd.DataFrame] = []  # To append df in loop
        _atom: dict[str, list[str]]  # Temporary dict to correct form
        symb_list: list[str] = []  # To save the symbol of each file
        style_list: list[str] = []  # To save the style of each forcefiled
        mix_list: list[str] = []  # To save the mix style of each forcefiled
        for f in self.param['files']:
            for atom in f['atoms']:
                _atom = {k: [v] for k, v in atom.items()}
                df = pd.DataFrame.from_dict(_atom, orient='columns')
                df_list.append(df)
                symb_list.append(f['symb'])
                style_list.append(f['style'])
                mix_list.append(f['mix'])
                del df
        param_df = pd.concat(df_list)
        param_df['f_symb'] = symb_list
        param_df['f_style'] = style_list
        param_df['mix'] = mix_list
        param_df.sort_values(by='type', axis=0, inplace=True)
        param_df.reset_index(inplace=True)
        param_df.drop(['index'], inplace=True, axis=1)
        param_df.index += 1
        return param_df

    def get_pairs(self) -> None:
        """find out the pairs that have somthing together"""
        self.param = self.mk_parameters()
        PARAMFIEL = 'parameters.lmp'
        print(f'{self.__class__.__name__}:\n'
              f'\tWritting: "{PARAMFIEL}"\n')
        self.__header: str  # str to write on top of each sections
        self.__header = ''.join(['#']*79)

        with open(PARAMFIEL, 'w') as f:
            try:
                pair_list = self.mk_pairs()
                self.write_pairs(pair_list, f)
            except KeyError:
                pass
            try:
                bond_df = self.mk_bond_couple()
                self.write_bond_copule(bond_df, f)
            except KeyError:
                pass
            try:
                angle_df = self.mk_angle_triple()
                self.write_angle_triple(angle_df, f)
            except KeyError:
                pass
            try:
                dihedral_df = self.mk_dihedral_quadruple()
                self.write_dihedral_qudrauple(dihedral_df, f)
            except KeyError:
                pass

    def write_pairs(self,
                    pair_list: list[tuple[int, int]],
                    f: typing.TextIO) -> None:
        """write pair nonbonding interactions"""
        _df = self.obj.Masses_df.copy()
        _df.index += 1
        # for pair in pair_list:
        style_set: set[str] = set(self.param['f_style'])
        styles = ' '.join(style_set)
        f.write(f'{self.__header}\n')
        f.write(f"\n")
        f.write(f'# interactions between pairs\n')
        f.write(f'\n')
        f.write(f'pair_style hybrid {styles}\n')
        f.write(f'\n')
        for i, pair in enumerate(pair_list):
            # print(i, self.param.iloc[pair[0]-1]['f_style'])
            name_i = _df['name'][pair[0]]
            name_j = _df['name'][pair[1]]
            file_i = self.param.iloc[pair[0]-1]['f_symb']
            file_j = self.param.iloc[pair[1]-1]['f_symb']
            if file_i is file_j:
                pair_style = self.param.iloc[pair[0]-1]['f_style']
                pair_style = self.drop_digit(pair_style)
            else:
                pair_style = 'lj/cut'
            args = []
            if pair[0] == pair[1]:
                sigma = self.param.iloc[pair[0]-1]['sigma']
                epsilon = self.param.iloc[pair[0]-1]['epsilon']
                r_cut = self.param.iloc[pair[0]-1]['r_cut']
            else:
                sigma = 's_mix'
                epsilon = 'e_mix'
                r_cut = 'r_mix'
            args.append(str(epsilon))
            args.append(str(sigma))
            args.append(str(r_cut))
            f.write(f'pair_coeff {pair[0]} {pair[1]}'
                    f' {pair_style} {" ".join(args)}')
            f.write(f' # {i+1} pair: {name_i} - {name_j}\n')
        f.write(f"\n")
        del _df

    def mix_geometric(self,
                      epsilon0: float,
                      epsilon1: float,
                      sigma0: float,
                      sigma1: float) -> tuple[float, float]:
        """mix interaction by geometric method form LAMMMPS manual"""
        epsilon = np.sqrt(epsilon0*epsilon1)
        sigma = np.sqrt(sigma0*sigma1)
        return epsilon, sigma
    
    def mix_arithmetic(self,
                      epsilon0: float,
                      epsilon1: float,
                      sigma0: float,
                      sigma1: float) -> tuple[float, float]:
        """mix interaction by arithmetic method form LAMMMPS manual"""
        epsilon = np.sqrt(epsilon0*epsilon1)
        sigma = 0.5*(sigma0+sigma1)
        return epsilon, sigma
    
    def mix_sixthpower(self,
                      epsilon0: float,
                      epsilon1: float,
                      sigma0: float,
                      sigma1: float) -> tuple[float, float]:
        """mix interaction by sixthpower method form LAMMMPS manual"""
        epsilon = 2 * np.sqrt(epsilon0*epsilon1)*sigma0**3*sigma1**3
        epsilon /= (sigma0**6 +sigma1**6)
        sigma = (0.5*(sigma0**6 + sigma**6))**(1/6)
        return epsilon, sigma


    def drop_digit(self, s: str) -> str:
        """drop numbers from string"""
        return re.sub(r" \d+", " ", s).strip()

    def write_bond_copule(self,
                          df: pd.DataFrame,
                          f: typing.TextIO) -> None:
        """write the bond pair and coefficents of the force field"""
        f.write(f"{self.__header}\n")
        f.write(f"\n")
        f.write(f"# coefficents for bonds interactions\n")
        f.write(f"\n")
        f.write(f"bond_style hybrid [args...]\n")
        f.write(f"\n")
        for i in range(len(df)):
            f.write(f'bond_coeff {df.iloc[i]["bond_typ"]}'
                    f' [style] [args]'
                    f' # bond_coeff for {df.iloc[i]["ai_name"]} -'
                    f' {df.iloc[i]["aj_name"]}')
            f.write(f"\n")
        f.write(f"\n")

    def write_angle_triple(self,
                           df: pd.DataFrame,
                           f: typing.TextIO) -> None:
        """write the triple atoms interaction that shares angle"""
        f.write(f"{self.__header}\n")
        f.write(f"\n")
        f.write(f"# coefficents for angles interactions\n")
        f.write(f"\n")
        f.write(f"angle_style hybrid [args...]\n")
        f.write(f"\n")
        for i in range(len(df)):
            f.write(f'angle_coeff {df.iloc[i]["angle_typ"]}'
                    f' [style] [args]'
                    f' # angle_coeff for {df.iloc[i]["ai_name"]} -'
                    f' {df.iloc[i]["aj_name"]} -'
                    f' {df.iloc[i]["ak_name"]}')
            f.write(f"\n")
        f.write(f"\n")

    def write_dihedral_qudrauple(self,
                                 df: pd.DataFrame,
                                 f: typing.TextIO) -> None:
        """write the quadrauple atoms interaction that shares dihedrals"""
        f.write(f"{self.__header}\n")
        f.write(f"\n")
        f.write(f"# coefficents for dihedrals interactions\n")
        f.write(f"\n")
        f.write(f"dihedral_style hybrid [args...]\n")
        f.write(f"\n")
        for i in range(len(df)):
            f.write(f'dihedral_coeff {df.iloc[i]["dihedral_typ"]}'
                    f' [style] [args]'
                    f' # dihedral_coeff for {df.iloc[i]["ai_name"]} -'
                    f' {df.iloc[i]["aj_name"]} -'
                    f' {df.iloc[i]["ak_name"]} -'
                    f' {df.iloc[i]["ah_name"]}')
            f.write(f"\n")
        f.write(f"\n")

    def mk_pairs(self) -> list[tuple[int, int]]:
        # Make pair of all atom type
        type_list = self.obj.Masses_df['typ']
        pair_list = itertools.combinations_with_replacement(type_list, 2)
        return pair_list

    def mk_bond_couple(self) -> pd.DataFrame:
        """make a DataFrame for type and name of bond interactions"""
        # Clairfy bond interactions
        _df: pd.DataFrame  # Temperary Dataframe to group Bonds_df
        bond_list: list[int]  # List of bond types
        ai_type: list[int] = []  # List of atom types
        aj_type: list[int] = []  # List of atom types
        ai_name: list[str] = []  # List of atoms name
        aj_name: list[str] = []  # List of atoms name
        i_id: int  # To temp save the atoms type
        j_id: int  # To temp save the atoms type
        _df_bpair: pd.DataFrame  # To return bond pair information
        _df = self.obj.Bonds_df.copy()
        _df = _df.groupby(by='typ').min()
        _df = _df.reset_index()  # To make "typ" as a column
        bond_list = _df['typ']
        _df_bpair = pd.DataFrame(bond_list)
        _df_bpair = _df_bpair.rename(columns={'typ': 'bond_typ'})
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
        del _df, ai_type, aj_type, ai_name, aj_name
        return _df_bpair

    def mk_angle_triple(self) -> pd.DataFrame:
        """write angle interaction between atoms"""
        # Clairfy bond interactions
        _df: pd.DataFrame  # Temperary Dataframe to group Angels_df
        angle_list: list[int]  # List of angle types
        ai_type: list[int] = []  # List of atom types
        aj_type: list[int] = []  # List of atom types
        ak_type: list[int] = []  # List of atom types
        ai_name: list[str] = []  # List of atoms name
        aj_name: list[str] = []  # List of atoms name
        ak_name: list[str] = []  # List of atoms name
        i_id: int  # To temp save the atoms type
        j_id: int  # To temp save the atoms type
        k_id: int  # To temp save the atoms type
        _df_apair: pd.DataFrame  # To return angle information

        _df = self.obj.Angles_df.copy()
        _df = _df.groupby(by=['typ']).min()
        _df = _df.reset_index()  # To make "typ" as a column
        angle_list = _df['typ']
        _df_apair = pd.DataFrame(angle_list)
        _df_apair = _df_apair.rename(columns={'typ': 'angle_typ'})
        for ai, aj, ak in zip(_df['ai'], _df['aj'], _df['ak']):
            i_id = self.obj.Atoms_df.iloc[ai-1]['typ']
            j_id = self.obj.Atoms_df.iloc[aj-1]['typ']
            k_id = self.obj.Atoms_df.iloc[ak-1]['typ']
            ai_type.append(i_id)
            aj_type.append(j_id)
            ak_type.append(k_id)
            ai_name.append(
                self.obj.Masses_df.iloc[i_id-1]['name']
                )
            aj_name.append(
                self.obj.Masses_df.iloc[j_id-1]['name']
                )
            ak_name.append(
                self.obj.Masses_df.iloc[k_id-1]['name']
                )
        _df_apair['ai_typ'] = ai_type
        _df_apair['aj_typ'] = aj_type
        _df_apair['ak_typ'] = ak_type
        _df_apair['ai_name'] = ai_name
        _df_apair['aj_name'] = aj_name
        _df_apair['ak_name'] = ak_name
        del _df, ai_type, aj_type, ak_type, ai_name, aj_name, ak_name
        return _df_apair

    def mk_dihedral_quadruple(self) -> pd.DataFrame:
        """write dihedral interaction between atoms"""
        # Clairfy bond interactions
        _df: pd.DataFrame  # Temperary Dataframe to group Dihedrals_df
        dihedral_list: list[int]  # List of dihedral types
        ai_type: list[int] = []  # List of atom types
        aj_type: list[int] = []  # List of atom types
        ak_type: list[int] = []  # List of atom types
        ah_type: list[int] = []  # List of atom types
        ai_name: list[str] = []  # List of atoms name
        aj_name: list[str] = []  # List of atoms name
        ak_name: list[str] = []  # List of atoms name
        ah_name: list[str] = []  # List of atoms name
        i_id: int  # To temp save the atoms type
        j_id: int  # To temp save the atoms type
        k_id: int  # To temp save the atoms type
        h_id: int  # To temp save the atoms type
        _df_dpair: pd.DataFrame  # To return dihedral information
        try:
            _df = self.obj.Dihedrals_df.copy()
            _df = _df.groupby(by=['typ']).min()
            _df = _df.reset_index()  # To make "typ" as a column
            dihedral_list = _df['typ']
            _df_dpair = pd.DataFrame(dihedral_list)
            _df_dpair = _df_dpair.rename(columns={'typ': 'dihedral_typ'})
            for ai, aj, ak, ah in zip(_df['ai'], _df['aj'],
                                      _df['ak'], _df['ah']):
                i_id = self.obj.Atoms_df.iloc[ai-1]['typ']
                j_id = self.obj.Atoms_df.iloc[aj-1]['typ']
                k_id = self.obj.Atoms_df.iloc[ak-1]['typ']
                h_id = self.obj.Atoms_df.iloc[ah-1]['typ']
                ai_type.append(i_id)
                aj_type.append(j_id)
                ak_type.append(k_id)
                ah_type.append(h_id)
                ai_name.append(
                    self.obj.Masses_df.iloc[i_id-1]['name']
                    )
                aj_name.append(
                    self.obj.Masses_df.iloc[j_id-1]['name']
                    )
                ak_name.append(
                    self.obj.Masses_df.iloc[k_id-1]['name']
                    )
                ah_name.append(
                    self.obj.Masses_df.iloc[h_id-1]['name']
                    )
            _df_dpair['ai_typ'] = ai_type
            _df_dpair['aj_typ'] = aj_type
            _df_dpair['ak_typ'] = ak_type
            _df_dpair['ah_typ'] = ah_type
            _df_dpair['ai_name'] = ai_name
            _df_dpair['aj_name'] = aj_name
            _df_dpair['ak_name'] = ak_name
            _df_dpair['ah_name'] = ah_name
            del _df, ai_type, aj_type, ak_type, ah_type
            del ai_name, aj_name, ak_name, ah_name
            return _df_dpair
        except AttributeError:
            empty_list: list[typing.Any] = []
            _df_dpair = pd.DataFrame(empty_list)
            return _df_dpair

    def keep_letter(self, s: str) -> str:
        """keep letter and remove the rest"""
        return re.sub(r'[^a-zA-Z]', '', s)
