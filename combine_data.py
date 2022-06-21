import re
import sys
import pandas as pd
import read_lmp_data as mlmp  # My lammps
import typing


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
        self.NATomTyp = self.update_atom_typ()
        self.Nmols = self.update_atom_mol()
        self.max_z = self.recenter_atoms()
        self.update_atom_name()
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
        NATomTyp: int = 0
        for i, f in enumerate(self.f_list):
            if i == 0:
                NATomTyp = self.l_headers[f].NATomTyp
            elif i > 0 and i < len(self.f_list):
                self.l_atoms[f]['typ'] += NATomTyp
                NATomTyp += self.l_headers[f].NATomTyp
            if i+1 > len(self.f_list):
                break
        return NATomTyp

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

    def recenter_atoms(self) -> float:
        """make sure all the atoms are not overlapping after stacking"""
        max_z: float = 0
        for i, f in enumerate(self.f_list):
            self.l_atoms[f] = self.return_to_zero(self.l_atoms[f])
            if i == 0:
                max_z = self.l_atoms[f]['z'].max()
            elif i > 0 and i < len(self.f_list):
                self.l_atoms[f]['z'] += max_z
                max_z = self.l_atoms[f]['z'].max()
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
                Natoms = self.l_headers[f].NATOMS
            elif i > 0 and i < len(self.f_list):
                # Update the the index of the second file
                self.l_bonds[f]['ai'] += Natoms
                self.l_bonds[f]['aj'] += Natoms
                # Add the number of the atoms of the current file
                Natoms += self.l_headers[f].NATOMS
            if i+1 > len(self.f_list):
                break

    def update_bond_typ(self) -> int:
        """update the number of each type in atoms card"""
        NBondTyp: int = 0
        for i, f in enumerate(self.f_list):
            if i == 0:
                NBondTyp = self.l_headers[f].NBondTyp
            elif i > 0 and i < len(self.f_list):
                self.l_bonds[f]['typ'] += NBondTyp
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
        """make the bond DataFrame"""
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

    def set_df_lists(self) -> None:
        """Set lists to append DataFrame in it"""
        self.l_atoms: dict[str, pd.DataFrame] = dict()
        self.l_bonds: dict[str, pd.DataFrame] = dict()
        self.l_angles: dict[str, pd.DataFrame] = dict()
        self.l_dihedrals: dict[str, pd.DataFrame] = dict()
        self.l_headers: dict[str, mlmp.Header] = dict()

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
        bond_name: list[str] = [
            f"{self.Atoms.iloc[ai-1]['name']}-"
            f"{self.Atoms.iloc[aj-1]['name']}"
            for ai, aj in zip(_Bond['ai'], _Bond['aj'])]
        _Bond['cmt'] = ["#" for _ in _Bond.index]
        _Bond['bond'] = bond_name
        self.Bonds = _Bond
        self.NBonds = bonds.NBonds
        self.NBondType = bonds.NBondTypes
        del _Bond

    def set_angles(self) -> None:
        """get angles with updating it with names"""
        angles = Angles(self.l_angles, self.l_headers, self.f_list)
        angles.mk_angles()
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
        self.NAngles = angles.NAgnles
        self.NAngleType = angles.NAngleType
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
            f.write(f"Data file from {INFILE} for silica slab\n")
            f.write(f"\n")
            self.write_numbers(f)
            self.write_box(f)
            # self.write_masses(f)
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
        f.write(f"{self.xlim[0]:.3f} {self.xlim[1]:.3f} xlo xhi\n")
        f.write(f"{self.ylim[0]:.3f} {self.ylim[1]:.3f} ylo yhi\n")
        f.write(f"{self.zlim[0]:.3f} {self.zlim[1]:.3f} zlo zhi\n")
        f.write(f"\n")

    def write_masses(self, f: typing.TextIO) -> None:
        f.write(f"Masses\n")
        f.write(f"\n")
        for item in self.typ.iterrows():
            f.write(f"{int(item[1]['typ']): 3}"
                    f"  {item[1]['mass']:.5f}"
                    f"  # {item[0]}\n"
                    )
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
        columns = ['typ', 'ai', 'aj', 'ak', 'ak', 'cmt', 'dihedral']
        self.system.Dihedrals.to_csv(f, sep=' ', index=True, columns=columns,
                                     header=None)
        f.write(f"\n")


INFILE = sys.argv[1:]
system = Combine(INFILE)
system.mk_lmp_df()
LMPFILE = 'interface.data'
lmp = WriteLmp(system)
lmp.mk_lmp()
