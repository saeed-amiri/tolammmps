import typing
import pandas as pd
import numpy as np
from pprint import pprint

"""
Reading data file from Ole Nickel, received on Jun 01, 2022
[ checking PEP8
    ~/.local/bin/pycodestyle lmp_to_pdb.py
and typing"
    ~/.local/bin/mypy lmp_to_pdb.py
]
Data contains many atoms,
The atoms with types 1 to 4 are the ones that belong to the SiO2
input:
    - LAMMPS data file
output:
    - LAMMPS data file

"""


class FILEERROR:
    """
    there is problem in the header of the INFILE,
    maybe a long header!\n
    """


class HEADER:
    """
    read haeder data of the data file
    check the number of the lines, atom, bond ... informations
    get the box , pairs, ... coefficents
    Use this class to read the header of the file (LAMMPS data file),
    and the file should have Masses with their name specified after (#)
    e.g.:

        Masses

        1 1.008000 # H
        2 16.000000 # OH
        3 16.000000 # OB
        4 28.059999 # Si

    it will return a few attributes for the class if they existed:
    Masses, Pair, and Angel and Dihedral coefficients. And also the name
    of the atoms types.
    The class BODY needs' names' to read the data file.
    """

    def __init__(self) -> None:
        self.atomsLine: int = 0
        self.atomsLine = self.check_file()
        print(f'number of header lines: {self.atomsLine}\n')
        self.read_header()

    def check_file(self) -> int:
        """ Check header
        input:
            - INFILE (lammps data file)
        output:
            - number of header lines
        """
        # An integer to prevent over-reading in case of header bugs
        MAXHEADER: int = 1000
        # track the number of lines in the hedaer
        linecount: int = 0
        with open(INFILE, 'r') as f:
            while True:
                linecount += 1
                line = f.readline()
                if line.strip().startswith('Atoms'):
                    atomsLine = linecount
                    break
                if linecount > MAXHEADER:
                    err = FILEERROR()
                    exit(err.__doc__)
                if not line:
                    exit("wrong data file\n")
        return atomsLine

    def read_header(self) -> None:
        """read header to get all the available info
        Read header now and get the data
        """
        # Setting dictionaries to save data of each block in the header
        self.set_attrs()
        # Setting flags to save data correctly
        Masses, PairCoeff, BondCoeff, AngleCoeff, DihedralCoeff, Atoms\
            = False, False, False, False, False, False
        # Track the number of lines
        linecount: int = 0
        with open(INFILE, 'r') as f:
            while True:
                linecount += 1
                if linecount > self.atomsLine:
                    break
                line: str = f.readline()
                if line.strip().endswith("atoms"):
                    self.NATOMS = int(line.strip().split(' ')[0])
                elif line.strip().endswith("atom types"):
                    self.NATomTyp = int(line.strip().split(' ')[0])
                elif line.strip().endswith("bonds"):
                    self.NBonds = int(line.strip().split(' ')[0])
                elif line.strip().endswith("bond types"):
                    self.NBondTyp = int(line.strip().split(' ')[0])
                elif line.strip().endswith("angles"):
                    self.NAngles = int(line.strip().split(' ')[0])
                elif line.strip().endswith("angle types"):
                    self.NAngleTyp = int(line.strip().split(' ')[0])
                elif line.strip().endswith("dihedrals"):
                    self.NDihedrals = int(line.strip().split(' ')[0])
                elif line.strip().endswith("dihedral typss"):
                    self.NDihedralsTyp = int(line.strip().split(' ')[0])
                elif line.strip().endswith("xhi"):
                    self.Xlim = self.get_axis_lim(line.strip().split('xlo')[0])
                elif line.strip().endswith("yhi"):
                    self.Ylim = self.get_axis_lim(line.strip().split('ylo')[0])
                elif line.strip().endswith("zhi"):
                    self.Zlim = self.get_axis_lim(line.strip().split('zlo')[0])
                # setting up Flages for reading the cards of data in the file
                elif line.strip().startswith("Masses"):
                    Masses, PairCoeff, BondCoeff, AngleCoeff, DihedralCoeff,\
                        Atoms = True, False, False, False, False, False
                elif line.strip().startswith("Pair"):
                    Masses, PairCoeff, BondCoeff, AngleCoeff, DihedralCoeff,\
                        Atoms = False, True, False, False, False, False
                elif line.strip().startswith("Bond Coeffs"):
                    Masses, PairCoeff, BondCoeff, AngleCoeff, DihedralCoeff,\
                        Atoms = False, False, True, False, False, False
                elif line.strip().startswith("Angle Coeffs"):
                    Masses, PairCoeff, BondCoeff, AngleCoeff, DihedralCoeff,\
                        Atoms = False, False, False, True, False, False
                elif line.strip().startswith("Dihedral Coeffs"):
                    Masses, PairCoeff, BondCoeff, AngleCoeff, DihedralCoeff,\
                        Atoms = False, False, False, False, True, False
                elif line.strip().startswith("Atoms"):
                    Masses, PairCoeff, BondCoeff, AngleCoeff, DihedralCoeff,\
                        Atoms = False, False, False, False, False, True
                elif line.strip():
                    if Masses:
                        self.get_masses(line.strip(), 'Masses')
                    elif PairCoeff:
                        self.get_pair_coeff(line.strip(), 'Pair')
                    elif BondCoeff:
                        self.get_bond_coeff(line.strip(), 'Bond')
                    elif AngleCoeff:
                        self.get_angle_coeff(line.strip(), 'Angle')
                    elif DihedralCoeff:
                        self.get_dihedral_coeff(line.strip(), 'Dihedral')
                if Atoms:
                    break
                if not line:
                    break

    def set_attrs(self) -> None:
        self.Names: dict[int, str] = dict()
        self.Masses: dict[int, float] = dict()
        self.PairCoeff: dict[int, typing.Any] = dict()
        self.BondCoeff: dict[int, typing.Any] = dict()
        self.AngleCoeff: dict[int, typing.Any] = dict()
        self.DihedralCoeff: dict[int, typing.Any] = dict()

    def get_axis_lim(self, lim) -> list:
        lim = lim.split(' ')
        lim = [float(item) for item in lim if item]
        return lim

    def get_masses(self, line, check) -> None:
        # stting the nth row of the dictionary
        if check not in line:
            typ = int(line.split(' ')[0])
            mass = float(line.split(' ')[1])
            atom_name = line.split('#')[1].strip()
            self.Masses[typ] = mass
            self.Names[typ] = atom_name

    def get_pair_coeff(self, line, check) -> None:
        # stting the nth row of the dictionary
        if check not in line:
            line = line.split(' ')
            typ = int(line[0])
            i_style = line[1]
            i_coeff = line[2:]
            self.PairCoeff[typ] = dict(
                                        style=i_style,
                                        coeff=i_coeff
                                       )

    def get_bond_coeff(self, line, check) -> None:
        # stting the nth row of the dictionary
        if check not in line:
            line = line.split(' ')
            typ = int(line[0])
            i_style = line[1]
            i_coeff = line[2:]
            self.BondCoeff[typ] = dict(
                                        style=i_style,
                                        coeff=i_coeff
                                       )

    def get_angle_coeff(self, line, check) -> None:
        # stting the nth row of the dictionary
        if check not in line:
            line = line.split(' ')
            typ = int(line[0])
            i_style = line[1]
            i_coeff = line[2:]
            self.AngleCoeff[typ] = dict(
                                        style=i_style,
                                        coeff=i_coeff
                                       )

    def get_dihedral_coeff(self, line, check) -> None:
        # stting the nth row of the dictionary
        if check not in line:
            line = line.split(' ')
            typ = int(line[0])
            i_style = line[1]
            i_coeff = line[2:]
            self.DihedralCoeff[typ] = dict(
                                            style=i_style,
                                            coeff=i_coeff
                                           )


class BODY:
    """
    read the data for atoms,velocities, bonds, angles, dihedrals
    It needs the names of the atoms read by HEADER class
    """

    def __init__(self, names) -> None:
        self.Name = names
        del names

    def read_body(self):
        self.Atoms, self.Velocities, self.Bonds, self.Angles, self.Dihedrals\
            = dict(), dict(), dict(), dict(), dict()
        Atoms, Velocities, Bonds, Angles, Dihedrals\
            = False, False, False, False, False

        with open(INFILE, 'r') as f:
            while True:
                line = f.readline()
                if line.strip().startswith('Atoms'):
                    Atoms, Velocities, Bonds, Angles, Dihedrals\
                        = True, False, False, False, False
                if line.strip().startswith('Velocities'):
                    Atoms, Velocities, Bonds, Angles, Dihedrals\
                        = False, True, False, False, False
                if line.strip().startswith('Bonds'):
                    Atoms, Velocities, Bonds, Angles, Dihedrals\
                        = False, False, True, False, False
                if line.strip().startswith('Angles'):
                    Atoms, Velocities, Bonds, Angles, Dihedrals\
                        = False, False, False, True, False
                if line.strip().startswith('Dihedrals'):
                    Atoms, Velocities, Bonds, Angles, Dihedrals\
                        = False, False, False, False, True
                if line.strip():
                    if Atoms:
                        self.get_atoms(line.strip())
                    if Velocities:
                        self.get_velocities(line.strip())
                    if Bonds:
                        self.get_bonds(line.strip())
                    if Angles:
                        self.get_angles(line.strip())
                    if Dihedrals:
                        self.get_dihedrals(line.strip())
                if not line:
                    break
            self.Atoms_df = pd.DataFrame.from_dict(
                            self.Atoms, orient='columns').T
            self.Bonds_df = pd.DataFrame.from_dict(self.Bonds).T

    def get_atoms(self, line) -> None:
        # stting the nth row of the dictionary
        if 'Atoms' not in line:
            line = line.split()
            line = [item for item in line if item]
            atom_id = int(line[0])
            i_mol = int(line[1])
            i_typ = int(line[2])
            i_charge = float(line[3])
            i_x = float(line[4])
            i_y = float(line[5])
            i_z = float(line[6])
            i_name = self.Name[i_typ]
            try:
                i_nx = int(line[7])
                i_ny = int(line[8])
                i_nz = int(line[9])
            except ValueError:
                i_nx = 0
                i_ny = 0
                i_nz = 0
            self.Atoms[atom_id] = dict(
                                        atom_id=atom_id,
                                        mol=i_mol,
                                        typ=i_typ,
                                        charge=i_charge,
                                        x=i_x,
                                        y=i_y,
                                        z=i_z,
                                        nx=i_nx,
                                        ny=i_ny,
                                        nz=i_nz,
                                        cmt='#',
                                        name=i_name
                                       )

    def get_velocities(self, line) -> None:
        # stting the nth row of the dictionary
        if 'Velocities' not in line:
            line = line.split()
            line = [item for item in line if item]
            atom_id = int(line[0])
            i_vx = float(line[1])
            i_vy = float(line[2])
            i_vz = float(line[3])
            self.Velocities[atom_id] = dict(
                                             vx=i_vx,
                                             vy=i_vy,
                                             vz=i_vz
                                            )

    def get_bonds(self, line) -> None:
        # stting the nth row of the dictionary
        if 'Bonds' not in line:
            line = line.split()
            line = [int(item) for item in line if item]
            bond_id = line[0]
            i_typ = int(line[1])
            i_ai = int(line[2])
            i_aj = int(line[3])
            self.Bonds[bond_id] = dict(
                                        typ=i_typ,
                                        ai=i_ai,
                                        aj=i_aj
                                       )

    def get_angles(self, line) -> None:
        # stting the nth row of the dictionary
        if "Angles" not in line:
            line = line.split()
            line = [int(item) for item in line if item]
            angle_id = line[0]
            i_typ = int(line[1])
            i_ai = int(line[2])
            i_aj = int(line[3])
            i_ak = int(line[4])
            self.Angles[angle_id] = dict(
                                         typ=i_typ,
                                         ai=i_ai,
                                         aj=i_aj,
                                         ak=i_ak
                                        )

    def get_dihedrals(self, line) -> None:
        # stting the nth row of the dictionary
        if "Dihedrals" not in line:
            line = line.split()
            line = [int(item) for item in line if item]
            dihedrals_id = line[0]
            i_typ = int(line[1])
            i_ai = int(line[2])
            i_aj = int(line[3])
            i_ak = int(line[4])
            i_ah = int(line[5])
            self.Dihedrals[dihedrals_id] = dict(
                                                 typ=i_typ,
                                                 ai=i_ai,
                                                 aj=i_aj,
                                                 ak=i_ak,
                                                 ah=i_ah
                                                )


class GeteSlab:
    """Get the infos for SiO2 slab only
    Slice the salb based fraction of the main data
    """

    def __init__(self, header, atoms) -> None:
        self.atoms = atoms
        self.header = header
        del atoms, header

    def get_slab(self) -> None:
        # Get the atoms coords based on the type
        self.get_atoms()
        # Get the bonds between all the SiO2 atoms
        self.get_bonds()
        self.get_infos()
        # Set min(x, y, z) -> (0, 0, 0), set the maximums
        self.move_to_center()

    def get_atoms(self) -> None:
        """extract eh SiO2 atoms and bonds from data file"""
        # Since the SiO2 atoms are seted from 1 to 4:
        atoms_df = self.atoms.Atoms_df[self.atoms.Atoms_df['typ'] < 5]
        self.atoms_df = self.crct_name(atoms_df)

    def get_bonds(self) -> None:
        bonds = []
        for i in range(self.header.NBonds):
            if self.atoms.Bonds_df.iloc[i]['ai'] in self.atoms_df['atom_id']:
                if self.atoms.Bonds_df.iloc[i][
                                        'aj'] in self.atoms_df['atom_id']:
                    bonds.append([self.atoms.Bonds_df.iloc[i]['typ'],
                                  atoms.Bonds_df.iloc[i]['ai'],
                                  atoms.Bonds_df.iloc[i]['aj']])
        self.bonds_df = pd.DataFrame(bonds, columns=['typ', 'ai', 'aj'])
        self.bonds_df.index += 1

    def move_to_center(self) -> None:
        # move to mins to orgin
        for i in ['x', 'y', 'z']:
            self.atoms_df[i] = self.atoms_df[i] - self.atoms_df[i].min()
        self.XMAX = np.max(self.atoms_df['x'])
        self.YMAX = np.max(self.atoms_df['y'])
        self.ZMAX = np.max(self.atoms_df['z'])

    def get_infos(self) -> None:
        self.get_name()
        self.count_atom()

    def get_name(self) -> None:
        # Get all the uniqe names
        self.Name = dict()
        for typ, name in self.atoms.Name.items():
            if typ < 5:
                self.Name[typ] = name

    def count_atom(self) -> None:
        # set the dictionary
        self.Count = dict()
        for typ, name in self.Name.items():
            self.Count[typ] = self.atoms_df.groupby(
                                ['name']).count()['x'].loc[name]

    def crct_name(self, df) -> pd.DataFrame:
        atom_name = []
        for i in range(len(df)):
            typ = df.iloc[i]['typ']
            q = df.iloc[i]['charge']
            if typ == 1:
                i_name = 'H'
            elif typ == 2:
                if q == -1.0:
                    i_name = 'OD'
                elif q == -0.9:
                    i_name = 'OND'
                elif q == -0.8:
                    i_name = 'OH'
                else:
                    i_name = None
            elif typ == 3:
                if q == -0.9:
                    i_name = 'OM'
                elif q == -0.8:
                    i_name = 'OB'
                else:
                    i_name = None
            elif typ == 4:
                if q == 1.5:
                    i_name = 'SM'
                elif q == 1.6:
                    i_name = 'Si'
            else:
                i_name = None
            atom_name.append(i_name)
        df['atom_name'] = atom_name
        return df


class CutSlab:
    """Cutting the slab
    The data from Ole is for 4*4 cell
    I need only one, here I ahve to cut the slab for it...
    A few things has to be consider:
        - Total charges of sliced sample
        - Ratio of particles compare to the orgiinal file

    !! The indicces hould be updated also the index for in the bonds
    also need to be updated !!
    """

    def __init__(self, atoms, bonds) -> None:
        self.atoms = atoms
        self.bonds = bonds
        self.x_cut_ratio = 0
        self.y_cut_ratio = 0
        del atoms, bonds

    def cut_slab(self) -> None:
        self.slice_slab(30.1, 30.42)
        self.set_new_index()
        self.set_index_dict()
        bonds = self.update_bond_index()
        self.bonds = self.mk_bonds_df(bonds)

    def slice_slab(self, xlim, ylim) -> None:
        if xlim:
            df_x = self.atoms[self.atoms['x'] < xlim]
        else:
            df_x = self.atoms
        if ylim:
            self.atoms = df_x[df_x['y'] < ylim]
        else:
            self.atoms = df_x

    def set_new_index(self) -> None:
        # update the index of atoms
        self.atoms = self.atoms.assign(index=pd.RangeIndex(start=1,
                                       stop=len(self.atoms)+1, step=1))
        self.atoms = self.atoms.set_index('index')

    def set_index_dict(self) -> None:
        """seting a dictionary to track the indexs and aupdate the
        bonds for atoms index
        """
        self.index_trace = {k: v for k, v in zip(self.atoms.atom_id,
                            self.atoms.index)}

    def update_bond_index(self) -> list:
        """Updating the indices in the bond
        and assign new index for them based on the self.index_trace"""
        updated_bonds = []
        new_bond_i = 0
        for i, (ai, aj) in enumerate(zip(self.bonds.ai, self.bonds.aj)):
            if aj in self.index_trace:
                if ai in self.index_trace:
                    try:
                        i_ai = self.index_trace[self.bonds.iloc[i]['ai']]
                        i_aj = self.index_trace[self.bonds.iloc[i]['aj']]
                        new_bond_i += 1
                        updated_bonds.append(
                            [new_bond_i, self.bonds.iloc[i]['typ'], i_aj,
                             i_ai])
                    except KeyError:
                        print(i, ai, aj, 'KeyError')
                    except IndexError:
                        print(i, ai, aj, 'IndexError')
        return updated_bonds

    def mk_bonds_df(self, bonds: list) -> pd.DataFrame:
        columns = ['id', 'typ', 'ai', 'aj']
        return pd.DataFrame(bonds, columns=columns)


class WriteData:
    """Write put LAMMPS data file
    write slab (whole or sliced) in LAMMPS version
    Masses section must be included
    INPUT:
        atoms and bonds from GetSlab class
    OUTPUT:
        LAMMPS data file
    """
    def __init__(self, atoms, bonds, Masses) -> None:
        self.atoms = atoms
        self.bonds = bonds
        self.Masses = Masses
        del atoms, bonds, Masses

    def write_lmp(self) -> None:
        """calling function to write data into a file"""
        # find box sizes
        self.set_box()
        # get number of atoms, types, bonds
        self.set_numbers()
        # write file
        self.write_data()
        print(self.atoms['charge'].sum())

    def set_box(self) -> None:
        """find Max and min of the data"""
        self.xlim = (self.atoms.x.min(), self.atoms.x.max())
        self.ylim = (self.atoms.y.min(), self.atoms.y.max())
        self.zlim = (self.atoms.z.min(), self.atoms.z.max())

    def set_numbers(self) -> None:
        """set the numbers of atoms, type, bonds"""
        self.Natoms = len(self.atoms)
        self.Nbonds = len(self.bonds)
        self.Natoms_type = np.max(self.atoms.typ)
        self.Nbonds_type = np.max(self.bonds.typ)

    def set_totals(self) -> None:
        """set the total numbers of charges, ..."""
        self.Tcharge = self.atoms['charge'].sum()

    def write_data(self) -> None:
        """write LAMMPS data file"""
        with open(OUTFILE, 'w') as f:
            f.write(f"Data file from Ole Nikle for silica slab\n")
            f.write(f"\n")
            self.write_numbers(f)
            self.write_box(f)
            self.write_masses(f)
            self.write_atoms(f)
            self.write_bonds(f)

    def write_numbers(self, f: typing.TextIO) -> None:
        f.write(f"{self.Natoms} atoms\n")
        f.write(f"{self.Natoms_type} atom types\n")
        f.write(f"{self.Nbonds} bonds\n")
        f.write(f"{self.Nbonds_type} bond types\n")
        f.write(f"\n")

    def write_box(self, f: typing.TextIO) -> None:
        f.write(f"{self.xlim[0]:.3f} {self.xlim[1]:.3f} xlo xhi\n")
        f.write(f"{self.ylim[0]:.3f} {self.ylim[1]:.3f} ylo yhi\n")
        f.write(f"{self.zlim[0]:.3f} {self.zlim[1]:.3f} zlo zhi\n")
        f.write(f"\n")

    def write_masses(self, f: typing.TextIO) -> None:
        f.write(f"Masses\n")
        f.write(f"\n")
        for k, v in self.Masses.items():
            if k < 5:
                f.write(f"{k} {v:.5f}\n")
        f.write(f"\n")

    def write_atoms(self, f: typing.TextIO) -> None:
        """Write atoms section inot file"""
        f.write(f"Atoms # full\n")
        f.write(f"\n")
        columns = ['mol', 'typ', 'charge', 'x', 'y', 'z', 'nx', 'ny', 'nz',
                   'cmt', 'atom_name']
        self.atoms.to_csv(f, sep=' ', index=True, columns=columns,
                          header=None)
        f.write(f"\n")
        f.write(f"\n")

    def write_bonds(self, f: typing.TextIO) -> None:
        f.write(f"Bonds\n")
        f.write(f"\n")
        columns = ['id', 'typ', 'ai', 'aj']
        self.bonds.to_csv(f, sep=' ', index=False, columns=columns,
                          header=None)


INFILE = 'merged.data'
OUTFILE = 'silica_ole_slic.data'
ole = HEADER()
atoms = BODY(ole.Names)
atoms.read_body()
slab = GeteSlab(ole, atoms)
slab.get_slab()
my_slice = CutSlab(slab.atoms_df, slab.bonds_df)
my_slice.cut_slab()
data = WriteData(my_slice.atoms, my_slice.bonds, ole.Masses)
data.write_lmp()
