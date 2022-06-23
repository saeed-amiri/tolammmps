import re
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
        - BonAnDi class instead of Bond, Angle, and Dihedral calsses,
        - Coniditoinal writing of the number of each cards,
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
        self.Natoms = self.update_atoms_df()
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
                Natoms = self.l_headers[f].NATOMS
            elif i > 0 and i < len(self.f_list):
                # Update the the index of the second file
                self.l_atoms[f]['atom_id'] += Natoms
                # Add the number of the atoms of the current file
                Natoms += self.l_headers[f].NATOMS
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

    def rm_special_str(self, char: list[typing.Any]) -> str:
        return re.sub('[^A-Za-z0-9]+', '', char)

    def stack_atoms(self) -> float:
        """make sure all the atoms are not overlapping after stacking"""
        max_z: float = 0
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
        """put all the staks in the center of mass of bottom one"""
        x_center: float = 0
        y_center: float = 0
        for i, f in enumerate(self.f_list):
            if i == 0:
                x_center += np.max(self.l_atoms[f]['x'])/2
                y_center += np.max(self.l_atoms[f]['y'])/2
            elif i > 0:
                _x_center: float = np.max(self.l_atoms[f]['x'])/2
                _y_center: float = np.max(self.l_atoms[f]['y'])/2
                self.l_atoms[f]['x'] += (x_center-_x_center)
                self.l_atoms[f]['y'] += (y_center-_y_center)


class Bonds:
    """update the bonds DataFrames and return uniq one"""

    def __init__(self,
                 l_bonds: dict[str, pd.DataFrame],
                 headers: dict[str, mlmp.Header],
                 f_list: list[str]) -> None:
        self.l_bonds = l_bonds
        self.l_headers = headers
        self.f_list = f_list
        del headers, l_bonds, f_list

    def mk_bonds(self) -> None:
        """make the bond DataFrame"""
        self.update_atoms_id()
        self.update_bond_typ()
        _Bonds: pd.DataFrame = self.append_bonds()
        _Bonds = _Bonds.reset_index()
        _Bonds.index += 1
        columns: list[str] = ['typ', 'ai', 'aj']
        self.Bonds: pd.DataFrame = _Bonds[columns].copy()
        self.NBonds: int = len(self.Bonds)
        self.NBondTypes: int = max(self.Bonds['typ'])
        del _Bonds

    def append_bonds(self) -> pd.DataFrame:
        """append bonds DataFrame with updataed id and type"""
        return pd.concat(self.l_bonds)

    def update_atoms_id(self) -> None:
        """Update the atom id (ai, aj)"""
        # Track the number of atom in each file
        Natoms: int = 0
        for i, f in enumerate(self.f_list):
            if i == 0:
                Natoms += self.l_headers[f].NATOMS
            elif i > 0 and i < len(self.f_list):
                try:
                    # Update the the index of the second file
                    self.l_bonds[f]['ai'] += Natoms
                    self.l_bonds[f]['aj'] += Natoms
                    # Add the number of the atoms of the current file
                    Natoms += self.l_headers[f].NATOMS
                except KeyError:
                    pass
            if i+1 > len(self.f_list):
                break

    def update_bond_typ(self) -> int:
        """update the number of each type in atoms card"""
        NBondTyp: int = 0
        for i, f in enumerate(self.f_list):
            if i == 0:
                NBondTyp = self.l_headers[f].NBondTyp
            elif i > 0 and i < len(self.f_list):
                try:
                    self.l_bonds[f]['typ'] += NBondTyp
                except KeyError:
                    pass
                NBondTyp += self.l_headers[f].NBondTyp
            if i+1 > len(self.f_list):
                break
        return NBondTyp


class Angles:
    """update the Angels DataFrames and return uniq one"""

    def __init__(self,
                 l_angles: dict[str, pd.DataFrame],
                 headers: dict[str, mlmp.Header],
                 f_list: list[str]) -> None:
        self.l_angles = l_angles
        self.l_headers = headers
        self.f_list = f_list
        del headers, l_angles, f_list

    def mk_angles(self) -> None:
        """make the angle DataFrame"""
        self.update_atoms_id()
        self.update_angle_typ()
        _Angles: pd.DataFrame = self.append_angles()
        _Angles = _Angles.astype(int)
        _Angles = _Angles.reset_index()
        _Angles.index += 1
        columns: list[str] = ['typ', 'ai', 'aj', 'ak']
        self.Angles: pd.DataFrame = _Angles[columns].copy()
        self.NAgnles: int = len(self.Angles)
        self.NAngleType: int = max(self.Angles['typ'])
        del _Angles

    def append_angles(self) -> pd.DataFrame:
        """append bonds DataFrame with updataed id and type"""
        return pd.concat(self.l_angles)

    def update_atoms_id(self) -> None:
        """Update the atom id (ai, aj, ak)"""
        # Track the number of atom in each file
        Natoms: int = 0
        for i, f in enumerate(self.f_list):
            if i == 0:
                Natoms = self.l_headers[f].NATOMS
            elif i > 0 and i < len(self.f_list):
                # Update the the index of the second file
                try:
                    self.l_angles[f]['ai'] += Natoms
                    self.l_angles[f]['aj'] += Natoms
                    self.l_angles[f]['ak'] += Natoms
                except KeyError:
                    pass
                # Add the number of the atoms of the current file
                Natoms += self.l_headers[f].NATOMS
            if i+1 > len(self.f_list):
                break

    def update_angle_typ(self) -> int:
        """update the number of each type in atoms card"""
        NAngleTyp: int = 0
        for i, f in enumerate(self.f_list):
            if i == 0:
                NAngleTyp = self.l_headers[f].NAngleTyp
            elif i > 0 and i < len(self.f_list):
                try:
                    self.l_angles[f]['typ'] += NAngleTyp
                except KeyError:
                    pass
                NAngleTyp += self.l_headers[f].NAngleTyp
            if i+1 > len(self.f_list):
                break
        return NAngleTyp


class Dihedrals:
    """update the Angels DataFrames and return uniq one"""
    def __init__(self,
                 l_dihedrals: dict[str, pd.DataFrame],
                 headers: dict[str, mlmp.Header],
                 f_list: list[str]) -> None:
        self.l_dihedrals = l_dihedrals
        self.l_headers = headers
        self.f_list = f_list
        del headers, l_dihedrals, f_list

    def mk_dihedrals(self) -> None:
        """make the bond DataFrame"""
        self.update_atoms_id()
        self.update_dihedral_typ()
        _Dihedrals: pd.DataFrame = self.append_dihedrals()
        _Dihedrals = _Dihedrals.astype(int)
        _Dihedrals = _Dihedrals.reset_index()
        _Dihedrals.index += 1
        columns: list[str] = ['typ', 'ai', 'aj', 'ak', 'ah']
        self.Dihedrals = _Dihedrals[columns]
        self.NDihedrals: int = len(self.Dihedrals['typ'])
        self.NDihedralType: int = max(self.Dihedrals['typ'])
        del _Dihedrals

    def append_dihedrals(self) -> pd.DataFrame:
        """append bonds DataFrame with updataed id and type"""
        return pd.concat(self.l_dihedrals)

    def update_atoms_id(self) -> None:
        """Update the atom id (ai, aj, ak, ah)"""
        # Track the number of atom in each file
        Natoms: int = 0
        for i, f in enumerate(self.f_list):
            if i == 0:
                Natoms = self.l_headers[f].NATOMS
            elif i > 0 and i < len(self.f_list):
                # Update the the index of the second file
                try:
                    self.l_dihedrals[f]['ai'] += Natoms
                    self.l_dihedrals[f]['aj'] += Natoms
                    self.l_dihedrals[f]['ak'] += Natoms
                    self.l_dihedrals[f]['ah'] += Natoms
                except KeyError:
                    pass
                # Add the number of the atoms of the current file
                Natoms += self.l_headers[f].NATOMS
            if i+1 > len(self.f_list):
                break

    def update_dihedral_typ(self) -> int:
        """update the number of each type in atoms card"""
        NDihedralTyp: int = 0
        for i, f in enumerate(self.f_list):
            if i == 0:
                NDihedralTyp = self.l_headers[f].NDihedralTyp
            elif i > 0 and i < len(self.f_list):
                try:
                    self.l_dihedrals[f]['typ'] += NDihedralTyp
                except KeyError:
                    pass
                NDihedralTyp += self.l_headers[f].NDihedralTyp
            if i+1 > len(self.f_list):
                break
        return NDihedralTyp


class Mass:
    """update the Mass card"""
    def __init__(self,
                 headers: dict[str, mlmp.Body],
                 f_list: list[str]) -> None:
        self.l_headers: dict[str, mlmp.Body] = headers
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

    def rm_special_str(self, char: list[typing.Any]) -> str:
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
        self.l_masses: mlmp.Header = dict()

    def get_data(self) -> None:
        """loop over all the files and make several DataFrame"""
        for f in self.f_list:
            print(f)
            o = mlmp.Header(f)
            self.l_headers[f] = o
            atoms = mlmp.Body(o.Names, f)
            atoms.read_body()
            self.l_atoms[f] = atoms.Atoms_df
            self.l_bonds[f] = atoms.Bonds_df
            self.l_angles[f] = atoms.Angles_df
            self.l_dihedrals[f] = atoms.Dihedrals_df
            self.l_masses[f] = o.Masses

    def set_atoms(self) -> None:
        """get atoms"""
        atoms = Atoms(self.l_atoms, self.l_headers, self.f_list)
        atoms.mk_atoms()
        self.Atoms = atoms.Atoms
        self.NAtoms: int = len(self.Atoms)
        self.NAtomType: int = max(self.Atoms['typ'])

    def set_bonds(self) -> None:
        """get bonds with updating it with names"""
        bonds = Bonds(self.l_bonds, self.l_headers, self.f_list)
        bonds.mk_bonds()
        _Bond: pd.DataFrame = bonds.Bonds
        _Bond = _Bond.astype(int)
        bond_name: list[str] = [
            f"{self.Atoms.iloc[ai-1]['name']}-"
            f"{self.Atoms.iloc[aj-1]['name']}"
            for ai, aj in zip(_Bond['ai'], _Bond['aj'])]
        _Bond['cmt'] = ["#" for _ in _Bond.index]
        _Bond['bond'] = bond_name
        self.Bonds = _Bond
        self.NBonds = int(bonds.NBonds)
        self.NBondType = int(bonds.NBondTypes)
        del _Bond

    def set_angles(self) -> None:
        """get angles with updating it with names"""
        angles = Angles(self.l_angles, self.l_headers, self.f_list)
        angles.mk_angles()
        self.NAngles = angles.NAgnles
        self.NAngleType = angles.NAngleType
        if self.NAngles > 0:
            _Angle: pd.DataFrame = angles.Angles
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
        dihedrals = Dihedrals(self.l_dihedrals, self.l_headers, self.f_list)
        dihedrals.mk_dihedrals()
        _Dihedrals: pd.DataFrame = dihedrals.Dihedrals
        angle_name: list[str] = [
            f"{self.Atoms.iloc[ai-1]['name']}-"
            f"{self.Atoms.iloc[aj-1]['name']}-"
            f"{self.Atoms.iloc[ak-1]['name']}-"
            f"{self.Atoms.iloc[ah-1]['name']}"
            for ai, aj, ak, ah in
            zip(_Dihedrals['ai'], _Dihedrals['aj'],
                _Dihedrals['ak'], _Dihedrals['ak'])]
        _Dihedrals['cmt'] = ["#" for _ in _Dihedrals.index]
        _Dihedrals['dihedral'] = angle_name
        self.Dihedrals = _Dihedrals
        self.NDihedrals = dihedrals.NDihedrals
        self.NDihedralType = dihedrals.NDihedralType
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
        self.write_data()
        # print(self.atoms['charge'].sum())

    def set_box(self) -> None:
        """find Max and min of the data"""
        self.xlim = (self.system.Atoms.x.min(), self.system.Atoms.x.max())
        self.ylim = (self.system.Atoms.y.min(), self.system.Atoms.y.max())
        self.zlim = (self.system.Atoms.z.min(), self.system.Atoms.z.max())

    def set_totals(self) -> None:
        """set the total numbers of charges, ..."""
        self.Tcharge = self.system.Atoms['charge'].sum()

    def write_data(self) -> None:
        """write LAMMPS data file"""
        with open(LMPFILE, 'w') as f:
            f.write(f"Data file from {INFILE} for an interface\n")
            f.write(f"\n")
            self.write_numbers(f)
            self.write_box(f)
            self.write_masses(f)
            self.write_atoms(f)
            self.write_bonds(f)
            self.write_angles(f)
            self.write_dihedrals(f)

    def write_numbers(self, f: typing.TextIO) -> None:
        f.write(f"{self.system.NAtoms} atoms\n")
        f.write(f"{self.system.NAtomType} atom types\n")
        f.write(f"{self.system.NBonds} bonds\n")
        f.write(f"{self.system.NBondType} bond types\n")
        f.write(f"{self.system.NAngles} angles\n")
        f.write(f"{self.system.NAngleType} angle types\n")
        f.write(f"{self.system.NDihedrals} dihedrals\n")
        f.write(f"{self.system.NDihedralType} dihedral types\n")
        f.write(f"\n")

    def write_box(self, f: typing.TextIO) -> None:
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
        f.write(f"Bonds\n")
        f.write(f"\n")
        columns = ['typ', 'ai', 'aj', 'cmt', 'bond']
        self.system.Bonds.to_csv(f, sep=' ', index=True, columns=columns,
                                 header=None)
        f.write(f"\n")

    def write_angles(self, f: typing.TextIO) -> None:
        f.write(f"Angles\n")
        f.write(f"\n")
        columns = ['typ', 'ai', 'aj', 'ak', 'cmt', 'angle']
        self.system.Angles.to_csv(f, sep=' ', index=True, columns=columns,
                                  header=None)
        f.write(f"\n")

    def write_dihedrals(self, f: typing.TextIO) -> None:
        f.write(f"Dihedrals\n")
        f.write(f"\n")
        columns = ['typ', 'ai', 'aj', 'ak', 'ah', 'cmt', 'dihedral']
        self.system.Dihedrals.to_csv(f, sep=' ', index=True, columns=columns,
                                     header=None)
        f.write(f"\n")


INFILE = sys.argv[1:]
system = Combine(INFILE)
system.mk_lmp_df()
LMPFILE = 'interface.data'
lmp = WriteLmp(system)
lmp.mk_lmp()
