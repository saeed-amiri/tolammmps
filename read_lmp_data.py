import sys
import typing
import pandas as pd


class FileErr:
    """
    there is problem in the header of the INFILE,
    maybe a long header!\n
    """


class Header:
    """
    read haeder of the data file
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

    def __init__(self, infile) -> None:
        self.infile: str = infile
        print(f'{self.__class__.__name__}:\n'
              f'\tReading: `{self.infile}`\n')
        self.atomsLine: int
        self.atomsLine = self.check_file()
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
        with open(self.infile, 'r') as f:
            while True:
                linecount += 1
                line = f.readline()
                if line.strip().startswith('Atoms'):
                    atomsLine = linecount
                    break
                if linecount > MAXHEADER:
                    err = FileErr()
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
        self.set_attr_zero()
        # Setting flags to save data correctly
        Masses, PairCoeff, BondCoeff, AngleCoeff, DihedralCoeff, Atoms\
            = False, False, False, False, False, False
        # Track the number of lines
        linecount: int = 0
        with open(self.infile, 'r') as f:
            while True:
                linecount += 1
                if linecount > self.atomsLine:
                    break
                line: str = f.readline()
                if line.strip().endswith("atoms"):
                    self.NAtoms = int(line.strip().split(' ')[0])

                elif line.strip().endswith("atom types"):
                    self.NAtomTyp = int(line.strip().split(' ')[0])

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

                elif line.strip().endswith("dihedral types"):
                    self.NDihedralTyp = int(line.strip().split(' ')[0])

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

    def set_attr_zero(self) -> None:
        """set the intial values to zero"""
        self.NAtoms = 0
        self.NBonds = 0
        self.NAngles = 0
        self.NAtomTyp = 0
        self.NBondTyp = 0
        self.NAngleTyp = 0
        self.NDihedrals = 0
        self.NDihedralTyp = 0

    def get_axis_lim(self, lim: typing.Any) -> list:
        lim = lim.split(' ')
        lim = [float(item) for item in lim if item]
        return lim

    def get_masses(self, line: str, check: str) -> None:
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


class Body(Header):
    """
    read the data for atoms, velocities, bonds, angles, dihedrals
    It needs the names of the atoms read by HEADER class
    """

    def __init__(self, infile) -> None:
        super().__init__(infile)
        self.infile: typing.IO = infile
        self.read_body()

    def read_body(self):
        self.Atoms, self.Velocities, self.Bonds, self.Angles, self.Dihedrals\
            = dict(), dict(), dict(), dict(), dict()
        Atoms, Velocities, Bonds, Angles, Dihedrals\
            = False, False, False, False, False

        with open(self.infile, 'r') as f:
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
            self.Angles_df = pd.DataFrame.from_dict(self.Angles).T
            self.Dihedrals_df = pd.DataFrame.from_dict(self.Dihedrals).T
            self.Masses_df = self.set_masses()
            del self.Atoms, self.Bonds, self.Angles, self.Dihedrals

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
            i_name = self.Names[i_typ]
            try:
                i_nx = int(line[7])
                i_ny = int(line[8])
                i_nz = int(line[9])
            except (ValueError, IndexError) as e:
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
        cmt_flag: bool = False
        if 'Bonds' not in line:
            if '#' in line:
                cmt_flag = True
            line = line.split()
            line = [item for item in line]
            line = [item for item in line if item]
            line[:4] = [int(item) for item in line[:4]]
            bond_id = line[0]
            i_typ = int(line[1])
            i_ai = int(line[2])
            i_aj = int(line[3])
            if cmt_flag:
                i_cmt = line[4]
                i_name = line[5]
                self.Bonds[bond_id] = dict(
                                           typ=i_typ,
                                           ai=i_ai,
                                           aj=i_aj,
                                           cmt=i_cmt,
                                           name=i_name
                                           )
            else:
                self.Bonds[bond_id] = dict(
                                           typ=i_typ,
                                           ai=i_ai,
                                           aj=i_aj
                                           )

    def get_angles(self, line) -> None:
        # stting the nth row of the dictionary
        if "Angles" not in line:
            if '#' in line:
                line = line.split('#')[0]
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
            if '#' in line:
                line = line.split('#')[0]
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

    def set_masses(self) -> pd.DataFrame:
        names_list: list[str] = []  # list to store all the names
        cmt_list: list[str] = []  # list to store '#'
        columns: list[str] = ['mass']  # column name of the DFs
        Masses_df = pd.DataFrame.from_dict(self.Masses,
                                           orient='index', columns=columns)
        Masses_df['typ'] = Masses_df.index
        for k, v in self.Masses.items():
            names_list.append(self.Names[k])
            cmt_list.append('#')
        Masses_df['cmt'] = cmt_list
        Masses_df['name'] = names_list
        return Masses_df


class ReadData(Body):
    """reading the input file
    This class call all other classes and make one output file
    """
    def __init__(self, infile) -> None:
        super().__init__(infile)


if __name__ == '__main__':
    ReadData(sys.argv[1])
