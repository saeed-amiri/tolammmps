import sys
import pandas as pd
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
        self.Atoms = self.append_atoms()
        print(self.Natoms, self.NATomTyp)

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
        self.Bonds = self.append_bonds()

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
        self.Angles = self.append_angles()

    def append_angles(self) -> pd.DataFrame:
        """append bonds DataFrame with updataed id and type"""
        return pd.concat(self.l_angles)

    def update_atoms_id(self) -> None:
        """Update the atom id (ai, aj)"""
        # Track the number of atom in each file
        Natoms: int = 0
        for i, f in enumerate(self.f_list):
                if i == 0:
                    Natoms = self.l_headers[f].NATOMS
                    print(f)
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
                try:
                    NAngleTyp = self.l_headers[f].NAngleTyp
                except KeyError:
                    pass
            elif i > 0 and i < len(self.f_list):
                try: 
                    self.l_angles[f]['typ'] += NAngleTyp
                    NAngleTyp += self.l_headers[f].NAngleTyp
                except KeyError:
                    pass
            if i+1 > len(self.f_list):
                break
        return NAngleTyp


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
        atoms = Atoms(self.l_atoms, self.l_headers, self.f_list)
        atoms.mk_atoms()
        self.Atoms = atoms.Atoms
        bonds = Bonds(self.l_bonds, self.l_headers, self.f_list)
        bonds.mk_bonds()
        self.Bonds = bonds.Bonds
        angles = Angles(self.l_angles, self.l_headers, self.f_list)
        angles.mk_angles()
        self.Angles = angles.Angles

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


INFILE = sys.argv[1:]
system = Combine(INFILE)
system.mk_lmp_df()
system.Atoms.to_csv('atoms', sep='\t', index=False, header=None)
system.Bonds.to_csv('bonds', sep='\t', index=False, header=None)
system.Angles.to_csv('angles', sep='\t', index=False, header=None)
