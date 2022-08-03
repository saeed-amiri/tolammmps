import sys
import typing
import numpy as np
import pandas as pd


class Doc:
    """write LAMMPS data file in 'full' atom style
    Input:
        class object with DataFrame of Atoms, Bonds, Angles, Dihedrals.
        Other number, will be find by this module itself.
    Output:
        LAMMPS data file
     """


class GetData:
    """Finding the number of each sections
    Input:
        DataFrames
    Output:
        Class object
    """
    def __init__(self, obj) -> None:
        self.get_atoms(obj.Atoms_df)
        if not obj.Bonds_df.empty:
            self.get_bonds(obj.Bonds_df)
        try:
            self.get_angles(obj.Angles_df)
        except AttributeError:
            pass
        try:
            self.get_dihedrals(obj.Dihedrals_df)
        except AttributeError:
            pass

    def get_atoms(self, df: pd.DataFrame) -> None:
        """get the number of atoms and thier types"""
        self.Natom_types: int  # Number of atom types
        self.Natoms: int  # Return len of df = Number of atoms
        self.xlo: float  # Low  limit of x column
        self.xhi: float  # High limit of x column
        self.ylo: float  # Low  limit of y column
        self.yhi: float  # High limit of y column
        self.zlo: float  # Low  limit of z column
        self.zhi: float  # High limit of z column

        self.Natom_types = max(df.typ)
        self.Natoms = len(df)
        self.xlo = np.min(df.x)
        self.xhi = np.max(df.x)
        self.ylo = np.min(df.y)
        self.yhi = np.max(df.y)
        self.zlo = np.min(df.z)
        self.zhi = np.max(df.z)

    def get_bonds(self, df: pd.DataFrame) -> None:
        """get the bond information"""
        self.Nbond_types: int  # Number of bond types
        self.Nbonds: int  # Return len of df = Number of bonds

        self.Nbond_types = np.max(df.typ)
        self.Nbonds = len(df)

    def get_angles(self, df: pd.DataFrame) -> None:
        """get the angle information"""
        self.Nangle_types: int  # Number of angle types
        self.Nangles: int  # Return len of df = Number of angles

        self.Nangle_types = np.max(df.typ)
        self.Nangles = len(df)

    def get_dihedrals(self, df: pd.DataFrame) -> None:
        """get the dihedrals information"""
        self.Ndihedral_types: int  # Number of dihedral types
        self.Ndihedrals: int  # Return len of df = Number of dihedrals

        self.Ndihedral_types = np.max(df.typ)
        self.Ndihedrals = len(df)


class WriteLmp(GetData):
    """write data in LAMMPS in full atom style"""
    def __init__(self, obj, output: str = 'blocked.data') -> None:
        super().__init__(obj)
        self.obj = obj
        self.fname = output
        print(f'{self.__class__.__name__}:\n'
              f'\tWrite `{self.fname}`\n')

    def write_lmp(self) -> None:
        """call all the function"""
        with open(self.fname, 'w') as f:
            self.write_header(f)
            self.write_body(f)

    def write_header(self, f: typing.TextIO) -> None:
        """write header of the file, including:
            Comments section: Names of the input files, date, code, dir
            Numbers section
            Mass section (if available)
            Box section
        """
        self.write_comments(f)
        self.write_numbers(f)
        self.write_box(f)
        self.write_masses(self.obj.Masses_df, f)

    def write_body(self, f: typing.TextIO) -> None:
        """write the body of the data file, including:
            atoms, bonds, angles, dihedrals
        """
        self.write_atoms(self.obj.Atoms_df, f)
        try:
            self.write_bonds(self.obj.Bonds_df, f)
        except AttributeError:
            pass
        try:
            self.write_angles(self.obj.Angles_df, f)
        except AttributeError:
            pass
        try:
            self.write_dihedrals(self.obj.Dihedrals_df, f)
        except AttributeError:
            pass

    def write_comments(self, f: typing.TextIO) -> None:
        """write comments on the top of the file"""
        f.write(f'# LAMMPS data file from: {sys.argv[1]} by {sys.argv[0]}\n')
        f.write(f'\n')

    def write_numbers(self, f: typing.TextIO) -> None:
        """write numbers of atoms, ..."""
        f.write(f'{self.Natoms} atoms\n')
        f.write(f'{self.Natom_types} atom types\n')
        try:
            f.write(f'{self.Nbonds} bonds\n')
            f.write(f'{self.Nbond_types} bond types\n')
        except AttributeError:
            pass
        try:
            f.write(f'{self.Nangles} angles\n')
            f.write(f'{self.Nangle_types} angle types\n')
        except AttributeError:
            pass
        try:
            f.write(f'{self.Ndihedral_types} dihedral types\n')
            f.write(f'{self.Ndihedrals} dihedrals\n')
        except AttributeError:
            pass
        f.write(f'\n')

    def write_masses(self, df: pd.DataFrame, f: typing.TextIO) -> None:
        """write mass section"""
        columns: list[str] = ['typ', 'mass', 'cmt', 'name']
        f.write(f"\n")
        f.write(f"Masses\n")
        f.write(f"\n")
        df.to_csv(f, sep=' ', index=False, columns=columns, header=None)
        f.write(f"\n")

    def write_box(self, f: typing.TextIO) -> None:
        """write box limits"""
        f.write(f'{self.xlo: 8.3f} {self.xhi+2: 8.3f} xlo xhi\n')
        f.write(f'{self.ylo: 8.3f} {self.yhi+2: 8.3f} ylo yhi\n')
        f.write(f'{self.zlo: 8.3f} {self.zhi+2: 8.3f} zlo zhi\n')
        f.write(f'\n')

    def write_atoms(self, df: pd.DataFrame, f: typing.TextIO) -> None:
        """write Atoms # full section"""
        if not df.empty:
            f.write(f'Atoms # full\n')
            f.write(f'\n')
            columns = ['atom_id', 'mol', 'typ', 'charge', 'x', 'y', 'z',
                       'nx', 'ny', 'nz', 'cmt', 'name']
            df = df.astype({'x': float, 'y':  float, 'z': float})
            df.to_csv(f, sep=' ', index=False, columns=columns, header=None,
                      float_format='%.8f')
            f.write(f'\n')
        else:
            exit(f'{self.__class__.__name__}\n'
                 f'\tError: Atoms section is empty\n')

    def write_bonds(self, df: pd.DataFrame, f: typing.TextIO) -> None:
        """write bonds section"""
        if not df.empty:
            f.write(f'Bonds\n')
            f.write(f'\n')
            columns = ['typ', 'ai', 'aj', 'cmt', 'name']
            df.to_csv(f, sep=' ', index=True, columns=columns, header=None)
            f.write(f'\n')
        else:
            print(f'{self.__class__.__name__}\n'
                  f'\tWARNING: Bonds section is empty\n')

    def write_angles(self, df: pd.DataFrame, f: typing.TextIO) -> None:
        """write angles section"""
        if not df.empty:
            f.write(f'Angles\n')
            f.write(f'\n')
            columns = ['typ', 'ai', 'aj', 'ak', 'cmt', 'name']
            df.to_csv(f, sep=' ', index=True, columns=columns, header=None)
            f.write(f'\n')
        else:
            print(f'{self.__class__.__name__}\n'
                  f'\tWARNING: Angels section is empty\n')

    def write_dihedrals(self, df: pd.DataFrame, f: typing.TextIO) -> None:
        """write dihedrals section"""
        if not df.empty:
            f.write(f'Dihedrals\n')
            f.write(f'\n')
            columns = ['typ', 'ai', 'aj', 'ak', 'ah', 'cmt', 'name']
            df.to_csv(f, sep=' ', index=True, columns=columns, header=None)
            f.write(f'\n')
        else:
            print(f'{self.__class__.__name__}\n'
                  f'\tWARNING: Dihedrals section is empty\n')
