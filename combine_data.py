import re
import sys
import typing
import itertools
import numpy as np
import pandas as pd
import read_lmp_data as mlmp  # My lammps reader


class Doc:
    """combining LAMMPS data file to prepare interface
    Input:
        LAMMPS data file[s]
    Output:
        LAMMPS data file
    [ checking PEP8
        ~/.local/bin/pycodestyle combine_data.py
        and typing"
        ~/.local/bin/mypy combine_data.py
    ]
    The order of data is based on the input order.
    For now, all the data only stack on top of each other in the
    z-direction.

    Usage:
        combine_data.py data1 data2 data3 ...

    !!!!!!!!!
    For now, all the input files should have different names and data.
    Files with the same names end up with wrong cards, mass, and types.
    If you use duplicate data files with different names, each file
    will get different types; Then, in the end, it has more complex
    pair interaction definitions.

    To fix:
        - BoAnDi class instead of Bond, Angle, and Dihedral calsses,
            Done Jun 23 2022
        - Coniditoinal writing of the number of each cards,
            Done Jun 24 2022
        - Writing LAMPPS as a modual to be used all the scripts,
        - The posibality of duplicate files,

    """


class Atoms:
    """ update and make a DataFrame for all atoms to be write into file
    All the inforamtion from atoms shuld come with this class:
    Natoms, box sizes, ...
    Input:
        A list of all the atom DataFrame and Headers
    """

    def __init__(self,
                 l_atoms: dict[str, pd.DataFrame],
                 headers: dict[str, mlmp.Header],
                 f_list: list[str]) -> None:
        self.l_atoms = l_atoms
        self.l_headers = headers
        self.f_list = f_list
        del headers, l_atoms, f_list

    def mk_atoms(self) -> None:
        """make atoms DataFrame"""
        self.NatomsNAtoms = self.update_atoms_df()
        self.NAtomTyp = self.update_atom_typ()
        self.Nmols = self.update_atom_mol()
        self.max_z = self.stack_atoms()
        self.update_atom_name()
        self.recenter_stak()
        Atoms = self.append_atoms()
        columns: list[str] = [
            'atom_id', 'mol', 'typ', 'charge', 'x', 'y', 'z',
            'nx', 'ny', 'nz', 'cmt', 'name']
        Atoms = Atoms.reset_index()
        self.Atoms = Atoms[columns].copy()

    def append_atoms(self) -> pd.DataFrame:
        """append atoms DataFrame with updataed id and type"""
        return pd.concat(self.l_atoms)

    def update_atoms_df(self) -> int:
        """Update the type and index of the atoms with number of next
        atoms in the next DataFrame
        """
        # Track the number of atom in each file
        Natoms: int = 0
        for i, f in enumerate(self.f_list):
            if i == 0:
                # jump the first file, no need for update
                # save the number of atoms
                Natoms = self.l_headers[f].NAtoms
            elif i > 0 and i < len(self.f_list):
                # Update the the index of the second file
                self.l_atoms[f]['atom_id'] += Natoms
                # Add the number of the atoms of the current file
                Natoms += self.l_headers[f].NAtoms
            if i+1 > len(self.f_list):
                break
        return Natoms

    def update_atom_typ(self) -> int:
        """update the number of each type in atoms card"""
        NAtomTyp: int = 0
        for i, f in enumerate(self.f_list):
            if i == 0:
                NAtomTyp = self.l_headers[f].NAtomTyp
            elif i > 0 and i < len(self.f_list):
                self.l_atoms[f]['typ'] += NAtomTyp
                NAtomTyp += self.l_headers[f].NAtomTyp
            if i+1 > len(self.f_list):
                break
        return NAtomTyp

    def update_atom_mol(self) -> int:
        """update the number of molecules"""
        Nmols: int = 0
        for i, f in enumerate(self.f_list):
            if i == 0:
                Nmols = self.l_atoms[f]['mol'].max()
            elif i > 0 and i < len(self.f_list):
                self.l_atoms[f]['mol'] += Nmols
                Nmols += self.l_atoms[f]['mol'].max()
            if i+1 > len(self.f_list):
                break
        return Nmols

    def update_atom_name(self) -> None:
        """Correct the name of the atoms"""
        for i, f in enumerate(self.f_list):
            # Drop the chars
            self.l_atoms[f]['name'] = [
                self.rm_special_str(item) for item in self.l_atoms[f]['name']]
            # Drop empty spaces
            self.l_atoms[f]['name'] = [
                item.strip() for item in self.l_atoms[f]['name']]
            # Keep only the first two letters
            self.l_atoms[f]['name'] = [
                item[0:2] for item in self.l_atoms[f]['name']]

    def rm_special_str(self, char: str) -> str:
        return re.sub('[^A-Za-z0-9]+', '', char)

    def stack_atoms(self) -> float:
        """make sure all the atoms are not overlapping after stacking"""
        max_z: float
        VACUME: int = 2  # Spcace between layers
        for i, f in enumerate(self.f_list):
            self.l_atoms[f] = self.return_to_zero(self.l_atoms[f])
            if i == 0:
                max_z = self.l_atoms[f]['z'].max()
                max_z += VACUME
            elif i > 0 and i < len(self.f_list):
                self.l_atoms[f]['z'] += max_z
                max_z = self.l_atoms[f]['z'].max()
                max_z += VACUME
            if i+1 > len(self.f_list):
                break
        return max_z

    def return_to_zero(self, df: pd.DataFrame) -> pd.DataFrame:
        """return minimums to zero in all three direction"""
        axes: list[str] = ['x', 'y', 'z']
        for ax in axes:
            min_x: float = df[ax].min()
            df[ax] -= min_x
        return df

    def recenter_stak(self) -> None:
        """put all the staks COM in the center of mass of bottom one"""
        x_center: float = 0  # X center of bottem layer
        y_center: float = 0  # Y center of bottem layer
        _x_center: float     # X center of other layers
        _y_center: float     # Y center of other layers
        for i, f in enumerate(self.f_list):
            if i == 0:
                x_center += np.max(self.l_atoms[f]['x'])/2
                y_center += np.max(self.l_atoms[f]['y'])/2
            elif i > 0:
                _x_center = np.max(self.l_atoms[f]['x'])/2
                _y_center = np.max(self.l_atoms[f]['y'])/2
                self.l_atoms[f]['x'] += (x_center-_x_center)
                self.l_atoms[f]['y'] += (y_center-_y_center)


class BoAnDi:
    """update the bonds, angles, and dihedrlas DataFrames and return
    unique one
    """

    def __init__(self,
                 l_df: dict[str, pd.DataFrame],
                 headers: dict[str, mlmp.Header],
                 f_list: list[str],
                 _type: str) -> None:
        self.l_df = l_df
        self.l_headers = headers
        self.f_list = f_list
        self._type = _type
        del headers, l_df, f_list, _type

    def mk_df(self) -> None:
        """make the bond DataFrame"""
        _columns: list[str] = self.set_columns()
        self.update_atoms_id(_columns)
        self.update_type()
        _df: pd.DataFrame = self.append_df()
        _df = _df.reset_index()
        _df.index += 1
        self.df: pd.DataFrame = _df[_columns].copy()
        self.Number: int = len(self.df)
        self.Ntype: int = max(self.df['typ'])
        del _df

    def set_columns(self) -> list[str]:
        """set a list for columns based on the _type"""
        _columns: list[str]
        if self._type == 'bond':
            _columns = ['typ', 'ai', 'aj']
        elif self._type == 'angle':
            _columns = ['typ', 'ai', 'aj', 'ak']
        elif self._type == 'dihedral':
            _columns = ['typ', 'ai', 'aj', 'ak', 'ah']
        return _columns

    def append_df(self) -> pd.DataFrame:
        """append bonds DataFrame with updataed id and type"""
        return pd.concat(self.l_df)

    def update_atoms_id(self, _columns: list[str]) -> None:
        """Update the atom id (ai, aj)"""
        # Track the number of atom in each file
        Natoms: int
        for i, f in enumerate(self.f_list):
            if i == 0:
                Natoms = self.l_headers[f].NAtoms
            elif i > 0 and i < len(self.f_list):
                try:
                    # Update the the index of the second file
                    for a in _columns:
                        if a != 'typ':
                            self.l_df[f][a] += Natoms
                    # Add the number of the atoms of the current file
                    Natoms += self.l_headers[f].NAtoms
                except KeyError:
                    pass
            if i+1 > len(self.f_list):
                break

    def update_type(self) -> int:
        """update the number of each type in atoms card"""
        Ntype: int   # To return
        _Ntype: int  # Temporary
        for i, f in enumerate(self.f_list):
            if self._type == 'bond':
                _Ntype = self.l_headers[f].NBondTyp
            elif self._type == 'angle':
                _Ntype = self.l_headers[f].NAngleTyp
            elif self._type == 'dihedral':
                _Ntype = self.l_headers[f].NDihedralTyp
            if i == 0:
                Ntype = _Ntype
            elif i > 0 and i < len(self.f_list):
                try:
                    self.l_df[f]['typ'] += _Ntype
                except KeyError:
                    pass
                Ntype += _Ntype
            if i+1 > len(self.f_list):
                break
            _Ntype = 0
        return Ntype


class Mass:
    """update the Mass card"""
    def __init__(self,
                 headers: dict[str, mlmp.Header],
                 f_list: list[str]) -> None:
        self.l_headers: dict[str, mlmp.Header] = headers
        self.f_list: list[str] = f_list
        del headers, f_list

    def mk_mass(self) -> None:
        """make the mass card"""
        Ntype: int = 0
        self.l_df: dict[str, typing.Any] = dict()
        self.l_name: dict[str, typing.Any] = dict()
        for f in self.f_list:
            df = self.mk_mass_df(f)
            df['typ'] += Ntype
            df_name = self.mk_names(f)
            df_name['typ'] += Ntype
            Ntype += len(self.l_headers[f].Masses)
            self.l_df[f] = df
            self.l_name[f] = df_name
        Names = pd.concat(self.l_name)
        self.Masses = pd.concat(self.l_df)
        self.Masses['cmt'] = ["#" for _ in self.Masses.typ]
        self.Masses[['name', 'f_name']] = Names[['name', 'f_name']]

    def rm_special_str(self, char: str) -> str:
        return re.sub('[^A-Za-z0-9]+', '', char)

    def mk_names(self, f: str) -> pd.DataFrame:
        """make names for the atoms"""
        columns = ['typ', 'name']
        df_name = pd.DataFrame(
            self.l_headers[f].Names.items(), columns=columns
            )
        df_name['name'] = [
            self.rm_special_str(item) for item in df_name['name']
            ]
        df_name['name'] = [item.strip() for item in df_name['name']]
        df_name['name'] = [item[0:2] for item in df_name['name']]
        df_name['f_name'] = [f"from_{f}" for _ in df_name['name']]
        return df_name

    def mk_mass_df(self, f: str) -> pd.DataFrame:
        """make mass dataframe for each file"""
        columns = ['typ', 'mass']
        df = pd.DataFrame(self.l_headers[f].Masses.items(), columns=columns)
        return df


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
        self.set_atoms()
        self.set_bonds()
        self.set_angles()
        self.set_dihedrals()
        self.set_masses()

    def set_df_lists(self) -> None:
        """Set lists to append DataFrame in it"""
        self.l_atoms: dict[str, pd.DataFrame] = dict()
        self.l_bonds: dict[str, pd.DataFrame] = dict()
        self.l_angles: dict[str, pd.DataFrame] = dict()
        self.l_dihedrals: dict[str, pd.DataFrame] = dict()
        self.l_headers: dict[str, mlmp.Header] = dict()
        self.l_masses: dict[str, dict[int, float]] = dict()

    def get_data(self) -> None:
        """loop over all the files and make several DataFrame"""
        for f in self.f_list:
            print(f)
            obj = mlmp.Header(f)
            self.l_headers[f] = obj
            atoms = mlmp.Body(obj.Names, f)
            atoms.read_body()
            self.l_atoms[f] = atoms.Atoms_df
            self.l_bonds[f] = atoms.Bonds_df
            self.l_angles[f] = atoms.Angles_df
            self.l_dihedrals[f] = atoms.Dihedrals_df
            self.l_masses[f] = obj.Masses

    def set_atoms(self) -> None:
        """get atoms"""
        atoms = Atoms(self.l_atoms, self.l_headers, self.f_list)
        atoms.mk_atoms()
        self.Atoms = atoms.Atoms
        self.NAtoms: int = len(self.Atoms)
        self.NAtomType: int = max(self.Atoms['typ'])

    def set_bonds(self) -> None:
        """get bonds with updating it with names"""
        bonds = BoAnDi(self.l_bonds, self.l_headers, self.f_list, 'bond')
        bonds.mk_df()
        bonds.df = bonds.df.astype(int)
        _Bond: pd.DataFrame = bonds.df
        self.NBonds = int(bonds.Number)
        self.NBondType = int(bonds.Ntype)
        if self.NBonds > 0:
            _Bond = _Bond.astype(int)
            bond_name: list[str] = [
                f"{self.Atoms.iloc[ai-1]['name']}-"
                f"{self.Atoms.iloc[aj-1]['name']}"
                for ai, aj in zip(_Bond['ai'], _Bond['aj'])]
            _Bond['cmt'] = ["#" for _ in _Bond.index]
            _Bond['bond'] = bond_name
            self.Bonds = _Bond
            del _Bond

    def set_angles(self) -> None:
        """get angles with updating it with names"""
        angles = BoAnDi(self.l_angles, self.l_headers, self.f_list, 'angle')
        angles.mk_df()
        angles.df = angles.df.astype(int)
        self.NAngles = int(angles.Number)
        self.NAngleType = int(angles.Ntype)
        if self.NAngles > 0:
            _Angle: pd.DataFrame = angles.df
            angle_name: list[str] = [
                f"{self.Atoms.iloc[ai-1]['name']}-"
                f"{self.Atoms.iloc[aj-1]['name']}-"
                f"{self.Atoms.iloc[ak-1]['name']}"
                for ai, aj, ak in
                zip(_Angle['ai'], _Angle['aj'], _Angle['ak'])]
            _Angle['cmt'] = ["#" for _ in _Angle.index]
            _Angle['angle'] = angle_name
            self.Angles = _Angle
            del _Angle

    def set_dihedrals(self) -> None:
        """get dihedrals with updating it with names"""
        dihedrals = BoAnDi(self.l_dihedrals, self.l_headers, self.f_list,
                           'dihedral')
        dihedrals.mk_df()
        dihedrals.df = dihedrals.df.astype(int)
        self.NDihedrals = int(dihedrals.Number)
        self.NDihedralType = int(dihedrals.Ntype)
        if self.NDihedrals > 0:
            _Dihedrals: pd.DataFrame = dihedrals.df
            dihedral_name: list[str] = [
                f"{self.Atoms.iloc[ai-1]['name']}-"
                f"{self.Atoms.iloc[aj-1]['name']}-"
                f"{self.Atoms.iloc[ak-1]['name']}-"
                f"{self.Atoms.iloc[ah-1]['name']}"
                for ai, aj, ak, ah in
                zip(_Dihedrals['ai'], _Dihedrals['aj'],
                    _Dihedrals['ak'], _Dihedrals['ak'])]
            _Dihedrals['cmt'] = ["#" for _ in _Dihedrals.index]
            _Dihedrals['dihedral'] = dihedral_name
            self.Dihedrals = _Dihedrals
            del _Dihedrals

    def set_masses(self) -> None:
        """make mass card"""
        mass = Mass(self.l_headers, self.f_list)
        mass.mk_mass()
        self.Masses = mass.Masses


class WriteLmp:
    """Write the data in a full atoms style for LAMMPS
    Input:
        atoms_df (DataFrame from PDBFILE: Pdb class)
        bonds_df, angles_df, dihedrals, impropers_df (DataFrame from
        PSFFILE Psf class)
    Output:
        A LAMMPS data file
    """

    def __init__(self, system: Combine) -> None:
        self.system = system
        print(f"Writting '{LMPFILE}' ...")
        del system

    def mk_lmp(self) -> None:
        """calling function to write data into a file"""
        # find box sizes
        self.set_box()
        # write file
        self.write_out()
        # print(self.atoms['charge'].sum())

    def set_box(self) -> None:
        """find Max and min of the data"""
        self.xlim = (self.system.Atoms.x.min(), self.system.Atoms.x.max())
        self.ylim = (self.system.Atoms.y.min(), self.system.Atoms.y.max())
        self.zlim = (self.system.Atoms.z.min(), self.system.Atoms.z.max())

    def set_totals(self) -> None:
        """set the total numbers of charges, ..."""
        self.Tcharge = self.system.Atoms['charge'].sum()

    def write_out(self) -> None:
        """Call write methods"""
        with open(LMPFILE, 'w') as f, open(PARAMFILE, 'w') as p:
            self.write_data(f)
            self.write_param(p)

    def write_data(self, f: typing.TextIO) -> None:
        """write LAMMPS data file"""
        f.write(f"Data file from {INFILE} for an interface\n")
        f.write(f"\n")
        self.write_numbers(f)
        self.write_box(f)
        self.write_masses(f)
        self.write_atoms(f)
        self.write_bonds(f)
        self.write_angles(f)
        self.write_dihedrals(f)

    def write_param(self, p: typing.TextIO) -> None:
        """write pair interactions file"""
        p.write(f"# Parameters file from {INFILE} for the ")
        p.write(f"interface file '{LMPFILE}'\n")
        p.write(f"\n")
        _header: list[str] = ['#']
        _header.extend(['*']*78)
        _header = ''.join(_header)
        # Making a DataFrame for bonds coeff
        _df = self.group_df(self.system.Bonds, 'bond')
        self.write_bond_couple(p, _df, _header)
        del _df
        # Making DataFrame for the pair the pair coeff
        _df = self.system.Masses.copy()
        _df = _df.reset_index()
        _df.index += 1
        self.write_atom_pair(p, _df, _header)
        del _df
        # Making a DataFrame for angle coeff
        _df = self.group_df(self.system.Angles, 'angle')
        self.write_angle_triple(p, _df, _header)
        # Making a DataFrame for dihedrals coeff
        _df = self.group_df(self.system.Dihedrals, 'dihedral')
        self.write_dihedrals_quadruple(p, _df, _header)
        del _df

    def group_df(self, df: pd.DataFrame, gb: str) -> pd.DataFrame:
        """Retrun grouped datafrma, gb: column name to group by"""
        _df = df.groupby(by=gb).min()
        _df = _df.reset_index()
        _df.index += 1
        _df = _df[gb]
        del df
        return _df

    def write_atom_pair(self,
                        p: typing.TextIO,
                        _df: pd.DataFrame,
                        _header: str) -> None:
        """Write pair interactions"""
        ll = list(_df['typ'])
        all_types = itertools.combinations_with_replacement(ll, 2)
        p.write(f"{_header}\n")
        p.write(f"# interactions between pairs\n")
        p.write(f"\n")
        p.write(f"pair_style hybrid [args...]\n")
        p.write(f"\n")
        for n, i in enumerate(all_types):
            name_i = _df['name'][i[0]]
            file_i = _df['f_name'][i[0]]
            name_j = _df['name'][i[1]]
            file_j = _df['f_name'][i[1]]
            p.write(f"#{n+1} pair: {name_i} {file_i} - {name_j} {file_j}\n")
            p.write(f"pair_coeff {i[0]} {i[1]} [pair_style] [args...] # ")
            p.write(f"{name_i} - {name_j}\n\n")
        p.write(f"\n")
        del _df

    def write_angle_triple(self,
                           p: typing.TextIO,
                           _df: pd.DataFrame,
                           _header: str) -> None:
        """Write angle interactions between partciles"""
        p.write(f"{_header}\n")
        p.write(f"# coefficents for angle interactions\n")
        p.write(f"\n")
        p.write(f"angle_style hybrid [args...]\n")
        p.write(f"\n")
        for n, i in enumerate(_df):
            p.write(f"#{n+1} angle_coeff for: {i}\n")
            p.write(f"angle_coeff {_df.index[n]} [style] [args...] \n")
            p.write(f"\n")
        p.write(f"\n")

    def write_dihedrals_quadruple(self,
                                  p: typing.TextIO,
                                  _df: pd.DataFrame,
                                  _header: str) -> None:
        """Write dihedral interactions"""
        p.write(f"{_header}\n")
        p.write(f"# coefficents for dihedral interactions\n")
        p.write(f"\n")
        p.write(f"dihedral_style hybrid [args...]\n")
        p.write(f"\n")
        for n, i in enumerate(_df):
            p.write(f"#{n+1} dihedral_coeff for: {i}\n")
            p.write(f"dihedral_coeff {_df.index[n]} [style] [args...] \n")
            p.write(f"\n")
        p.write(f"\n")

    def write_bond_couple(self,
                          p: typing.TextIO,
                          _df: pd.DataFrame,
                          _header: str) -> None:
        """Write coefficents for bonds"""
        p.write(f"{_header}\n")
        p.write(f"# coefficents for bonds interactions\n")
        p.write(f"\n")
        p.write(f"bond_style hybrid [args...]\n")
        p.write(f"\n")
        for n, i in enumerate(_df):
            p.write(f"#{n+1} bond_coeff for: {i}\n")
            p.write(f"bond_coeff {_df.index[n]} [style] [args...] \n")
            p.write(f"\n")
        p.write(f"\n")

    def write_numbers(self, f: typing.TextIO) -> None:
        """Write the number of each types"""
        f.write(f"{self.system.NAtoms} atoms\n")
        f.write(f"{self.system.NAtomType} atom types\n")
        if self.system.NBonds > 0:
            f.write(f"{self.system.NBonds} bonds\n")
            f.write(f"{self.system.NBondType} bond types\n")
        if self.system.NAngles > 0:
            f.write(f"{self.system.NAngles} angles\n")
            f.write(f"{self.system.NAngleType} angle types\n")
        if self.system.NDihedrals > 0:
            f.write(f"{self.system.NDihedrals} dihedrals\n")
            f.write(f"{self.system.NDihedralType} dihedral types\n")
        f.write(f"\n")

    def write_box(self, f: typing.TextIO) -> None:
        """Write the box sizes"""
        f.write(f"{self.xlim[0]-5:.3f} {self.xlim[1]+5:.3f} xlo xhi\n")
        f.write(f"{self.ylim[0]-5:.3f} {self.ylim[1]+5:.3f} ylo yhi\n")
        f.write(f"{self.zlim[0]-5:.3f} {self.zlim[1]+5:.3f} zlo zhi\n")
        f.write(f"\n")

    def write_masses(self, f: typing.TextIO) -> None:
        columns = ['typ', 'mass', 'cmt', 'name', 'f_name']
        f.write(f"Masses\n")
        f.write(f"\n")
        self.system.Masses.to_csv(f, sep=' ', index=False, columns=columns,
                                  header=None)
        f.write(f"\n")

    def write_atoms(self, f: typing.TextIO) -> None:
        """Write atoms section inot file"""
        f.write(f"Atoms # full\n")
        f.write(f"\n")
        columns = ['atom_id', 'mol', 'typ', 'charge', 'x', 'y', 'z', 'nx',
                   'ny', 'nz', 'cmt', 'name']
        self.system.Atoms = self.system.Atoms.astype(
            {'x': float, 'y':  float, 'z': float}
        )
        self.system.Atoms.to_csv(f, sep=' ', index=False, columns=columns,
                                 header=None, float_format='%.8f')
        f.write(f"\n")

    def write_bonds(self, f: typing.TextIO) -> None:
        if self.system.NBonds > 0:
            f.write(f"Bonds\n")
            f.write(f"\n")
            columns = ['typ', 'ai', 'aj', 'cmt', 'bond']
            self.system.Bonds.to_csv(f, sep=' ', index=True, columns=columns,
                                     header=None)
            f.write(f"\n")

    def write_angles(self, f: typing.TextIO) -> None:
        if self.system.NAngles > 0:
            f.write(f"Angles\n")
            f.write(f"\n")
            columns = ['typ', 'ai', 'aj', 'ak', 'cmt', 'angle']
            self.system.Angles.to_csv(f, sep=' ', index=True, columns=columns,
                                      header=None)
            f.write(f"\n")

    def write_dihedrals(self, f: typing.TextIO) -> None:
        if self.system.NDihedrals > 0:
            f.write(f"Dihedrals\n")
            f.write(f"\n")
            columns = ['typ', 'ai', 'aj', 'ak', 'ah', 'cmt', 'dihedral']
            self.system.Dihedrals.to_csv(f, sep=' ', index=True,
                                         columns=columns, header=None)
            f.write(f"\n")


INFILE = sys.argv[1:]
system = Combine(INFILE)
system.mk_lmp_df()
LMPFILE = 'interface.data'
PARAMFILE = 'parmeters.data'
lmp = WriteLmp(system)
lmp.mk_lmp()
