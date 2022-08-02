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


class MakeParamDf:
    """Make DataFrame of atoms, bonds, angles, dihedrals
    Input:
        obj from updated parameter-file
    Output:
        DataFrames
    """
    def __init__(self, param) -> None:
        self.param = param
        self.mk_dataframes()
        del param

    def mk_dataframes(self) -> None:
        """call all the methods"""
        self.lj_df = self.mk_lj_piar_df()
        self.mix_df = self.mk_mix_pair_df()  # Parameters for the i&j pairs
        self.bond_df = self.mk_bond_df()
        self.angle_df = self.mk_angle_df()
        self.dihedral_df = self.mk_dihedral_df()

    def mk_lj_piar_df(self) -> pd.DataFrame:
        """make dataframe from updated parameter file for pair intera-
        ctions.
        It is the main DataFrame read and updated from the input para-
        meter file.
        It contains all the information from the "atoms" section in
        the parameter input file.
        """
        df: pd.DataFrame  # Temporary dataframe
        lj_df: pd.DataFrame  # The main dataframe for pairwise interaction
        df_list: list[pd.DataFrame] = []  # To append df in loop
        _atom: dict[str, list[str]]  # Temporary dict to correct form
        symb_list: list[str] = []  # To save the symbol of each file
        style_list: list[str] = []  # To save the style of each forcefiled
        mix_list: list[str] = []  # To save the mix style of each forcefiled
        for f in self.param['files']:
            if 'atoms' in f:
                for atom in f['atoms']:
                    _atom = {k: [v] for k, v in atom.items()}
                    df = pd.DataFrame.from_dict(_atom, orient='columns')
                    df_list.append(df)
                    symb_list.append(f['symb'])
                    style_list.append(f['style'])
                    mix_list.append(f['mix'])
                    del df
            else:
                exit(f'\tERROR: there is no atom defeined in the file: '
                     f'`{f["file"]}`')
        lj_df = pd.concat(df_list)
        del df_list
        lj_df['f_symb'] = symb_list
        lj_df['f_style'] = style_list
        lj_df['mix'] = mix_list
        lj_df.sort_values(by='type', inplace=True, axis=0)
        lj_df.reset_index(inplace=True)
        lj_df.drop(['index'], inplace=True, axis=1)
        lj_df.index += 1
        return lj_df

    def mk_mix_pair_df(self) -> pd.DataFrame:
        """make a dataframe for mixed pair interactions
        It is DataFrame from the "Pairs" section of the parameter input
        file.
        It includes information on how the atoms in the different files
        interact.
        """
        df: pd.DataFrame  # Temporary dataframe
        pair_df: pd.DataFrame  # Return the of mixed pair interaction
        pair_dict: dict[str, list[str]]  # dict of the pairs
        pair_list: list[pd.DataFrame] = []  # list of the temporary df
        for p in self.param['Pairs']:
            pair_dict = {k: [v] for k, v in p.items()}
            df = pd.DataFrame.from_dict(pair_dict)
            pair_list.append(df)
        pair_df = pd.concat(pair_list)
        pair_df.reset_index(inplace=True)
        pair_df.drop(['index'], inplace=True, axis=1)
        pair_df.index += 1
        return pair_df

    def mk_bond_df(self) -> pd.DataFrame:
        """make dataframe for all the bonds information.
        It contains the bonds with updated "types"
        """
        df: pd.DataFrame  # Temporary dataframe
        bond_df: pd.DataFrame  # The main dataframe for bonds' info
        df_list: list[pd.DataFrame] = []  # To append df in loop
        _bond: dict[str, list[str]]  # Temporary dict to correct form
        symb_list: list[str] = []  # To save the symbol of each file
        for f in self.param['files']:
            for bond in f['bonds']:
                _bond = {k: [v] for k, v in bond.items()}
                df = pd.DataFrame.from_dict(_bond, orient='columns')
                df_list.append(df)
                symb_list.append(f['symb'])
                del df
        bond_df = pd.concat(df_list)
        del df_list
        bond_df['f_symb'] = symb_list
        bond_df.sort_values(by='type', inplace=True, axis=0)
        bond_df.reset_index(inplace=True)
        bond_df.drop(['index'], inplace=True, axis=1)
        bond_df.index += 1
        return bond_df

    def mk_angle_df(self) -> pd.DataFrame:
        """make dataframe for all the angles information.
        It contains the bonds with updated "types"
        """
        df: pd.DataFrame  # Temporary dataframe
        angle_df: pd.DataFrame  # The main dataframe for angles' info
        df_list: list[pd.DataFrame] = []  # To append df in loop
        _angle: dict[str, list[str]]  # Temporary dict to correct form
        symb_list: list[str] = []  # To save the symbol of each file
        for f in self.param['files']:
            if 'angles' in f:
                for angle in f['angles']:
                    _angle = {k: [v] for k, v in angle.items()}
                    df = pd.DataFrame.from_dict(_angle, orient='columns')
                    df_list.append(df)
                    symb_list.append(f['symb'])
                    del df
            else:
                pass
        angle_df = pd.concat(df_list)
        del df_list
        angle_df['f_symb'] = symb_list
        angle_df.sort_values(by='type', inplace=True, axis=0)
        angle_df.reset_index(inplace=True)
        angle_df.drop(['index'], inplace=True, axis=1)
        angle_df.index += 1
        return angle_df

    def mk_dihedral_df(self) -> pd.DataFrame:
        """make dataframe for all the dihedrals information.
        It contains the bonds with updated "types"
        """
        df: pd.DataFrame  # Temporary dataframe
        dihedral_df: pd.DataFrame  # The main dataframe for dihedrals' info
        df_list: list[pd.DataFrame] = []  # To append df in loop
        _dihedral: dict[str, list[str]]  # Temporary dict to correct form
        symb_list: list[str] = []  # To save the symbol of each file
        for f in self.param['files']:
            if 'dihedrals' in f:
                for dihedral in f['dihedrals']:
                    _dihedral = {k: [v] for k, v in dihedral.items()}
                    df = pd.DataFrame.from_dict(_dihedral, orient='columns')
                    df_list.append(df)
                    symb_list.append(f['symb'])
                    del df
            else:
                pass
        if df_list:
            dihedral_df = pd.concat(df_list)
            dihedral_df['f_symb'] = symb_list
            dihedral_df.sort_values(by='type', inplace=True, axis=0)
            dihedral_df.reset_index(inplace=True)
            dihedral_df.drop(['index'], inplace=True, axis=1)
            dihedral_df.index += 1
        else:
            dihedral_df = pd.DataFrame(df_list)
        del df_list
        return dihedral_df


class WriteParam(MakeParamDf):
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
        super().__init__(self.param)
        self.get_pairs()
        del obj, bs

    def get_pairs(self) -> None:
        """find out the pairs that have somthing together"""
        PARAMFIEL = 'parameters.lmp'
        print(f'{self.__class__.__name__}:\n'
              f'\tWritting: `{PARAMFIEL}`\n')
        self.__header: str  # str to write on top of each sections
        self.__header = ''.join(['#']*79)

        with open(PARAMFIEL, 'w') as f:
            try:
                pair_list = self.mk_pairs()
                self.write_pairs(pair_list, f)
            except KeyError:
                pass
            try:
                bond_pair_df = self.mk_bond_couple()
                self.write_bond_copule(bond_pair_df, f)
            except KeyError:
                pass
            try:
                angle_triple_df = self.mk_angle_triple()
                self.write_angle_triple(angle_triple_df, f)
            except KeyError:
                pass
            try:
                dihedral_quadruple_df = self.mk_dihedral_quadruple()
                self.write_dihedral_qudrauple(dihedral_quadruple_df, f)
            except KeyError:
                pass

    def write_pairs(self,
                    pair_list: list[tuple[int, int]],
                    f: typing.TextIO) -> None:
        """write pair nonbonding interactions"""
        _df = self.obj.Masses_df.copy()
        _df.index += 1
        style_set: set[str] = set(self.lj_df['f_style'])
        styles = ' '.join(style_set)
        f.write(f'{self.__header}\n')
        f.write(f"\n")
        f.write(f'# interactions between pairs\n')
        f.write(f'\n')
        f.write(f'pair_style hybrid {styles}\n')
        f.write(f'\n')
        for i, pair in enumerate(pair_list):
            name_i: str = _df['name'][pair[0]]
            name_j: str = _df['name'][pair[1]]
            file_i: str = self.lj_df.iloc[pair[0]-1]['f_symb']
            file_j: str = self.lj_df.iloc[pair[1]-1]['f_symb']
            epsilon_i: float = self.lj_df.iloc[pair[0]-1]['epsilon']
            epsilon_j: float = self.lj_df.iloc[pair[1]-1]['epsilon']
            sigma_i: float = self.lj_df.iloc[pair[0]-1]['sigma']
            sigma_j: float = self.lj_df.iloc[pair[1]-1]['sigma']
            r_cut_i: float = self.lj_df.iloc[pair[0]-1]['r_cut']
            r_cut_j: float = self.lj_df.iloc[pair[1]-1]['r_cut']
            mix_i: str = self.lj_df.iloc[pair[0]-1]['mix']
            mix_j: str = self.lj_df.iloc[pair[1]-1]['mix']

            pair_style: str = self.set_pair_style(file_i, file_j, pair)

            args: str = self.set_pair_args(
                epsilon_i, epsilon_j, sigma_i, sigma_j, r_cut_i, r_cut_j, pair,
                mix_i, mix_j, file_i, file_j)

            f.write(f'pair_coeff {pair[0]} {pair[1]}'
                    f' {pair_style} {args}')
            f.write(f'  # {i+1} pair: {name_i} - {name_j}\n')
        f.write(f"\n")
        del _df

    def set_pair_args(self,
                      epsilon_i: float,
                      epsilon_j: float,
                      sigma_i: float,
                      sigma_j: float,
                      r_cut_i: float,
                      r_cut_j: float,
                      pair: tuple[int, int],
                      mix_i: str,
                      mix_j: str,
                      file_i: str,
                      file_j: str) -> str:
        """make a sequnce of the interaction arguments"""
        mix: str = None  # "mix" style LAMMPS equation
        r_cut: float  # Calculated r_cut
        epsilon: float  # Calculated epsilon
        sigma: float  # Calculated sigma
        args: list[typing.Any] = []  # Make the arguments for interactions
        if pair[0] == pair[1]:
            sigma = sigma_i
            epsilon = epsilon_i
            r_cut = r_cut_i
        else:
            if file_i == file_j:  # if both atom are from the same file
                epsilon, sigma, r_cut = self.mixed_sigma_epsilon(
                    epsilon_i, epsilon_j, sigma_i, sigma_j,
                    r_cut_i, r_cut_j, mix_i)
            else:
                files = [file_i, file_j]
                check_mix = itertools.permutations(files, 2)
                for m in check_mix:
                    pair_mix = f"{m[0]}{m[1]}"
                    if pair_mix in list(self.mix_df['pair']):
                        mix = self.mix_df.loc[self.mix_df['pair'] ==
                                              pair_mix]['mix'][1]
                if mix is None:
                    exit(f'\tError!: The interactions between'
                         f' files are not defeiend')
                epsilon, sigma, r_cut = self.mixed_sigma_epsilon(
                    epsilon_i, epsilon_j, sigma_i, sigma_j,
                    r_cut_i, r_cut_j, mix)

        args.append(f'{epsilon: .5f}')
        args.append(f'{sigma: .5f}')
        # args.append(f'{r_cut: .5f}')
        return " ".join(args)

    def mixed_sigma_epsilon(self,
                            epsilon_i: float,
                            epsilon_j: float,
                            sigma_i: float,
                            sigma_j: float,
                            r_cut_i: float,
                            r_cut_j: float,
                            mix: str) -> tuple[float, float, float]:
        """return mixed sigma and epsilon for ij mixed pairs"""
        if mix == 'geometric':
            epsilon, sigma, r_cut = self.mix_geometric(
                epsilon_i, epsilon_j, sigma_i, sigma_j, r_cut_i, r_cut_j)
        elif mix == 'arithmetic':
            epsilon, sigma, r_cut = self.mix_arithmetic(
                epsilon_i, epsilon_j, sigma_i, sigma_j, r_cut_i, r_cut_j)
        elif mix == 'sixthpower':
            epsilon, sigma, r_cut = self.mix_sixthpower(
                epsilon_i, epsilon_j, sigma_i, sigma_j, r_cut_i, r_cut_j)
        return epsilon, sigma, r_cut

    def set_pair_style(self,
                       file_i: str,
                       file_j: str,
                       pair: tuple[int, int]) -> str:
        """set the pair style for each one of them"""
        if file_i is file_j:
            pair_style: str = self.lj_df.iloc[pair[0]-1]['f_style']
            pair_style = self.drop_digit(pair_style)
        else:
            pair_style = 'lj/cut'
        return pair_style

    def write_bond_copule(self,
                          df: pd.DataFrame,
                          f: typing.TextIO) -> None:
        """write the bond pair and coefficents of the force field"""
        bonds_set: set[str]  # Get the unique bond style in the files
        bonds_set = set(self.bond_df['style'])
        style_args: str  # A string contain all the bonds style
        style_args = " ".join(bonds_set)
        f.write(f'{self.__header}\n')
        f.write(f'\n')
        f.write(f'# coefficents for bonds interactions\n')
        f.write(f'\n')
        f.write(f'bond_style hybrid {style_args}\n')
        f.write(f'\n')
        for i in range(len(df)):
            style: str = self.bond_df.iloc[i]['style']
            args: str = self.bond_args(i)
            f.write(f'bond_coeff {df.iloc[i]["bond_typ"]}'
                    f' {style} {args}'
                    f' # bond_coeff for {df.iloc[i]["name"]}'
                    f' (name: {self.bond_df.iloc[i]["name"]})'
                    )
            f.write(f"\n")
        f.write(f"\n")

    def bond_args(self, i_loc: int) -> str:
        """return str contains arguments for the bond coeffs"""
        args: str
        args = f'{self.bond_df.iloc[i_loc]["kbond"]: .5f} ' \
               f'{self.bond_df.iloc[i_loc]["r"]: .5f}'
        return args

    def write_angle_triple(self,
                           df: pd.DataFrame,
                           f: typing.TextIO) -> None:
        """write the triple atoms interaction that shares angle"""
        angles_set: set[str]  # Get the unique angle style in the files
        angles_set = set(self.angle_df['style'])
        style_args: str  # A string contain all the angles style
        style_args = " ".join(angles_set)
        f.write(f'{self.__header}\n')
        f.write(f'\n')
        f.write(f'# coefficents for angles interactions\n')
        f.write(f'\n')
        f.write(f'angle_style hybrid {style_args}\n')
        f.write(f'\n')
        for i in range(len(df)):
            style: str = self.angle_df.iloc[i]['style']
            args: str = self.angle_args(i)
            f.write(f'angle_coeff {df.iloc[i]["angle_typ"]}'
                    f' {style} {args}'
                    f' # angle_coeff for {df.iloc[i]["name"]}'
                    f' (name: {self.angle_df.iloc[i]["name"]})'
                    )
            f.write(f"\n")
        f.write(f"\n")

    def angle_args(self, i_loc: int) -> str:
        """return str contains arguments for the angle coeffs"""
        args: str
        args = f'{self.angle_df.iloc[i_loc]["kangle"]: .5f} '\
               f'{self.angle_df.iloc[i_loc]["angle"]: .5f}'
        return args

    def write_dihedral_qudrauple(self,
                                 df: pd.DataFrame,
                                 f: typing.TextIO) -> None:
        """write the quadrauple atoms interaction that shares dihedrals"""
        dihedrals_set: set[str]  # Get the unique dihedral style in the files
        dihedrals_set = set(self.dihedral_df['style'])
        style_args: str  # A string contain all the dihedrals style
        style_args = " ".join(dihedrals_set)
        f.write(f"{self.__header}\n")
        f.write(f"\n")
        f.write(f"# coefficents for dihedrals interactions\n")
        f.write(f"\n")
        f.write(f"dihedral_style hybrid {style_args}\n")
        f.write(f"\n")
        for i in range(len(df)):
            style: str = self.dihedral_df.iloc[i]['style']
            args: str = self.dihedral_args(i)
            f.write(f'dihedral_coeff {df.iloc[i]["dihedral_typ"]}'
                    f' {style} {args}'
                    f' # dihedral_coeff for {df.iloc[i]["name"]}'
                    f' (name: {self.dihedral_df.iloc[i]["name"]})'
                    )
            f.write(f"\n")
        f.write(f"\n")

    def dihedral_args(self, i_loc: int) -> str:
        """return str contains arguments for the angle coeffs"""
        args: str
        args = f'{self.dihedral_df.iloc[i_loc]["k1"]: .5f} '\
               f'{self.dihedral_df.iloc[i_loc]["k2"]: .5f} '\
               f'{self.dihedral_df.iloc[i_loc]["k3"]: .5f} '\
               f'{self.dihedral_df.iloc[i_loc]["k4"]: .5f}'
        return args

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
        i_id: int  # To temp save the atoms type
        j_id: int  # To temp save the atoms type
        _df_bpair: pd.DataFrame  # To return bond pair information
        _df = self.obj.Bonds_df.copy()
        _df['ai'] = pd.to_numeric(_df['ai'])
        _df = _df.loc[_df.groupby(by=['typ'])['ai'].idxmin(), :]
        _df = _df.reset_index()  # To make "typ" as a column
        bond_list = _df['typ']
        _df_bpair = pd.DataFrame(bond_list)
        _df_bpair = _df_bpair.rename(columns={'typ': 'bond_typ'})
        _df_bpair['name'] = _df['name'].copy()
        for ai, aj in zip(_df['ai'], _df['aj']):
            i_id = self.obj.Atoms_df.iloc[ai-1]['typ']
            j_id = self.obj.Atoms_df.iloc[aj-1]['typ']
            ai_type.append(i_id)
            aj_type.append(j_id)
        _df_bpair['ai_typ'] = ai_type
        _df_bpair['aj_typ'] = aj_type
        del _df, ai_type, aj_type
        _df_bpair.index += 1
        return _df_bpair

    def mk_angle_triple(self) -> pd.DataFrame:
        """write angle interaction between atoms"""
        # Clairfy bond interactions
        _df: pd.DataFrame  # Temperary Dataframe to group Angels_df
        angle_list: list[int]  # List of angle types
        ai_type: list[int] = []  # List of atom types
        aj_type: list[int] = []  # List of atom types
        ak_type: list[int] = []  # List of atom types
        i_id: int  # To temp save the atoms type
        j_id: int  # To temp save the atoms type
        k_id: int  # To temp save the atoms type
        _df_apair: pd.DataFrame  # To return angle information
        _df = self.obj.Angles_df.copy()
        _df = _df.loc[_df.groupby(by=['typ'])['ai'].idxmin(), :]
        _df = _df.reset_index()  # To make "typ" as a column
        angle_list = _df['typ']
        _df_apair = pd.DataFrame(angle_list)
        _df_apair = _df_apair.rename(columns={'typ': 'angle_typ'})
        _df_apair['name'] = _df['name'].copy()
        for ai, aj, ak in zip(_df['ai'], _df['aj'], _df['ak']):
            i_id = self.obj.Atoms_df.iloc[ai-1]['typ']
            j_id = self.obj.Atoms_df.iloc[aj-1]['typ']
            k_id = self.obj.Atoms_df.iloc[ak-1]['typ']
            ai_type.append(i_id)
            aj_type.append(j_id)
            ak_type.append(k_id)
        _df_apair['ai_typ'] = ai_type
        _df_apair['aj_typ'] = aj_type
        _df_apair['ak_typ'] = ak_type
        del _df, ai_type, aj_type, ak_type
        _df_apair.index += 1
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
        i_id: int  # To temp save the atoms type
        j_id: int  # To temp save the atoms type
        k_id: int  # To temp save the atoms type
        h_id: int  # To temp save the atoms type
        _df_dpair: pd.DataFrame  # To return dihedral information
        try:
            _df = self.obj.Dihedrals_df.copy()
            _df = _df.loc[_df.groupby(by=['typ'])['ai'].idxmin(), :]
            _df = _df.reset_index()  # To make "typ" as a column
            dihedral_list = _df['typ']
            _df_dpair = pd.DataFrame(dihedral_list)
            _df_dpair = _df_dpair.rename(columns={'typ': 'dihedral_typ'})
            _df_dpair['name'] = _df['name'].copy()
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
            _df_dpair['ai_typ'] = ai_type
            _df_dpair['aj_typ'] = aj_type
            _df_dpair['ak_typ'] = ak_type
            _df_dpair['ah_typ'] = ah_type
            del _df, ai_type, aj_type, ak_type, ah_type
            del ai_name
            _df_dpair.index += 1
            return _df_dpair
        except AttributeError:
            empty_list: list[typing.Any] = []
            _df_dpair = pd.DataFrame(empty_list)
            _df_dpair.index += 1
            return _df_dpair

    def mix_geometric(self,
                      epsilon_i: float,
                      epsilon_j: float,
                      sigma_i: float,
                      sigma_j: float,
                      r_cut_i: float,
                      r_cut_j: float) -> tuple[float, float, float]:
        """mix interaction by geometric method form LAMMMPS manual"""
        epsilon = np.sqrt(epsilon_i*epsilon_j)
        sigma = np.sqrt(sigma_i*sigma_j)
        r_cut = np.sqrt(r_cut_i*r_cut_j)
        return epsilon, sigma, r_cut

    def mix_arithmetic(self,
                       epsilon_i: float,
                       epsilon_j: float,
                       sigma_i: float,
                       sigma_j: float,
                       r_cut_i: float,
                       r_cut_j: float) -> tuple[float, float, float]:
        """mix interaction by arithmetic method form LAMMMPS manual"""
        epsilon = np.sqrt(epsilon_i*epsilon_j)
        sigma = 0.5*(sigma_i+sigma_j)
        r_cut = 0.5*(r_cut_i+r_cut_j)
        return epsilon, sigma, r_cut

    def mix_sixthpower(self,
                       epsilon_i: float,
                       epsilon_j: float,
                       sigma_i: float,
                       sigma_j: float,
                       r_cut_i: float,
                       r_cut_j: float) -> tuple[float, float, float]:
        """mix interaction by sixthpower method form LAMMMPS manual"""
        epsilon: float = 2 * np.sqrt(epsilon_i*epsilon_j)*sigma_i**3*sigma_j**3
        epsilon /= (sigma_i**6 + sigma_j**6)
        sigma: float = (0.5*(sigma_i**6 + sigma_j**6))**(1/6)
        r_cut: float = (0.5*(r_cut_i**6 + r_cut_j**6))**(1/6)
        return epsilon, sigma, r_cut

    def drop_digit(self, s: str) -> str:
        """drop numbers from string"""
        # Drop dot '.' in case of a float number
        s = re.sub(r" *\.+", " ", s).strip()
        return re.sub(r" \d+", " ", s).strip()

    def keep_letter(self, s: str) -> str:
        """keep letter and remove the rest"""
        return re.sub(r'[^a-zA-Z]', '', s)
