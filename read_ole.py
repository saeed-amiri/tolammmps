import pandas as pd
from pprint import pprint
import time
import concurrent.futures

"""
Reading data file from Ole Nickel, received on Jun 01, 2022
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
    maybe long header!\n
    """


start = time.time()


class HEADER:
    """
    read haeder data of the data file
    check the number of the lines, atom, bond ... informations
    get the box , pairs, ... coefficents
    """

    def __init__(self) -> None:
        self.atomsLine = 0
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
        # An integer to prevent overreading in case of header bugs
        MAXHEADER = 1000
        # track the number of lines in the hedaer
        linecount = 0
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

    def read_header(self):
        """read header to get all the available info
        Read header now get data
        """
        # Setting dictionaries to save data of each block in the header
        self.Masses, self.PairCoeff, self.BondCoeff, self.AngleCoeff,\
            self.DihedralCoeff = dict(), dict(), dict(), dict(), dict()
        # Setting flags to save data correctly
        Masses, PairCoeff, BondCoeff, AngleCoeff, DihedralCoeff, Atoms\
            = False, False, False, False, False, False
        # Track the number of lines
        linecount = 0
        with open(INFILE, 'r') as f:
            while True:
                linecount += 1
                if linecount > self.atomsLine:
                    break
                line = f.readline()
                # self.process_type_box(line)
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
                elif line.strip().startswith("Masses"):
                    Masses, PairCoeff, BondCoeff, AngleCoeff, DihedralCoeff, Atoms\
                         = True, False, False, False, False, False
                elif line.strip().startswith("Pair"):
                    Masses, PairCoeff, BondCoeff, AngleCoeff, DihedralCoeff, Atoms\
                        = False, True, False, False, False, False
                elif line.strip().startswith("Bond Coeffs"):
                    Masses, PairCoeff, BondCoeff, AngleCoeff, DihedralCoeff, Atoms\
                        = False, False, True, False, False, False
                elif line.strip().startswith("Angle Coeffs"):
                    Masses, PairCoeff, BondCoeff, AngleCoeff, DihedralCoeff, Atoms\
                        = False, False, False, True, False, False
                elif line.strip().startswith("Dihedral Coeffs"):
                    Masses, PairCoeff, BondCoeff, AngleCoeff, DihedralCoeff, Atoms\
                        = False, False, False, False, True, False
                elif line.strip().startswith("Atoms"):
                    Masses, PairCoeff, BondCoeff, AngleCoeff, DihedralCoeff, Atoms\
                        = False, False, False, False, False, True
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

    def get_axis_lim(self, lim) -> list:
        lim = lim.split(' ')
        lim = [float(item) for item in lim if item]
        return lim

    def get_masses(self, line, check) -> dict:
        if check not in line:
            typ = line.split(' ')[0]
            mass = float(line.split(' ')[1])
            self.Masses[typ] = mass
        else:
            pass

    def get_pair_coeff(self, line, check)-> dict:
        if check not in line: 
            line = line.split(' ')
            typ = line[0]
            self.PairCoeff[typ]=dict(style=line[1], coeff=line[2:])
        else: pass
    
    def get_bond_coeff(self, line, check)-> dict:
        if check not in line:
            line = line.split(' ')
            typ = line[0]
            self.BondCoeff[typ]=dict(style=line[1], coeff=line[2:])
        else: pass
        
    def get_angle_coeff(self, line, check)-> dict:
        if check not in line:
            line = line.split(' ')
            typ = line[0]
            self.AngleCoeff[typ]=dict(style=line[1], coeff=line[2:])
        else: pass
    
    def get_dihedral_coeff(self, line, check)-> dict:
        if check not in line:
            line = line.split(' ')
            typ = line[0]
            self.DihedralCoeff[typ]=dict(style=line[1], coeff=line[2:])
        else: pass



class BODY:
    """
    read the data for atoms,velocities, bonds, angles, dihedrals
    """

    def __init__(self) -> None:
        pass
    
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
                    if Atoms :self.get_atoms(line.strip())
                    if Velocities: self.get_velocities(line.strip())
                    if Bonds: self.get_bonds(line.strip())
                    if Angles: self.get_angles(line.strip())
                    if Dihedrals: self.get_dihedrals(line.strip())
                if not line: break
            self.Atoms_df = pd.DataFrame.from_dict(
                            self.Atoms, orient='columns').T
            self.Atoms_df = self.Atoms_df.astype(
                            {'atom_id':'int','mol': 'int', 'typ': 'int'})
            self.Atoms_df = self.Atoms_df.astype(
                            {'nx':'int','ny': 'int', 'nz': 'int'})
            self.Bonds_df = pd.DataFrame.from_dict(self.Bonds).T
            self.Bonds_df = self.Bonds_df.astype({'ai': 'int', 'aj': 'int'})
    
    def get_atoms(self, line) -> dict:
        if 'Atoms' not in line:
            line = line.split(" ")
            line = [item for item in line if item]
            atom_id = int(line[0])
            i_mol=int(line[1])
            i_typ=int(line[2]); i_charg=float(line[3]) 
            i_x=float(line[4]);i_y=float(line[5]); i_z=float(line[6])
            try:
                i_nx=str(line[7]); i_ny=str(line[8]); i_nz=str(line[9])
            except:
                i_nx=0; i_ny=0; i_nz=0
                pass
            self.Atoms[atom_id]=dict(
                atom_id=atom_id ,mol=i_mol, typ=i_typ, charg=i_charg,\
                x=i_x,y=i_y, z=i_z, nx=i_nx, ny=i_ny, nz=i_nz)
        else: pass

    def get_velocities(self, line) -> dict:
        if 'Velocities' not in line:
            line = line.split(" ")
            line = [item for item in line if item]
            atom_id = int(line[0])
            i_vx = float(line[1]); i_vy = float(line[2]); i_vz = float(line[3])
            self.Velocities[atom_id] = dict(vx = i_vx, vy = i_vy, vz = i_vz)
        else: pass
    
    def get_bonds(self, line) -> dict:
        if 'Bonds' not in line:
            line = line.split(" ")
            line = [int(item) for item in line if item]
            bond_id = line[0]
            i_typ = line[1]; i_ai = line[2]; i_aj = line[3]
            self.Bonds[bond_id] = dict(typ = i_typ, ai = i_ai, aj = i_aj)
        else: pass

    def get_angles(self, line) -> dict:
        if "Angles" not in line:
            line = line.split(" ")
            line = [int(item) for item in line if item]
            angle_id = line[0]
            i_typ = line[1]; i_ai = line[2]; i_aj = line[3]; i_ak = line[4]
            self.Angles[angle_id] = dict(
                typ = i_typ, ai = i_ai, aj = i_aj, ak = i_ak)

    def get_dihedrals(self, line) -> dict:
        if "Dihedrals" not in line:
            line = line.split(" ")
            line = [int(item) for item in line if item]
            dihedrals_id = line[0]
            i_typ = line[1]
            i_ai = line[2]
            i_aj = line[3]
            i_ak = line[4]
            i_ah = line[5]
            self.Dihedrals[dihedrals_id] = dict(
                typ = i_typ, ai = i_ai, aj = i_aj, ak = i_ak, ah = i_ah)
INFILE = 'merged.data'

ole = HEADER()
atoms = BODY()
atoms.read_body()
MAXI = atoms.Atoms_df['x'].max()/4
MAXI = 25.7184*2
MAXI = 32.278
MAYI = atoms.Atoms_df['z'].max()/2
# MAYI = 22.5981
MAYI = 34.358

atoms_df0 = atoms.Atoms_df[atoms.Atoms_df['typ']<5 ]
atoms_df1 = atoms_df0[atoms_df0['x']<MAXI]
atoms_df = atoms_df1[atoms_df1['y']<MAYI]
bonds = []
start = time.time()
for i in range(ole.NBonds):
# def slice_df(i):
    if atoms.Bonds_df.iloc[i]['ai'] in atoms_df['atom_id']:
        if atoms.Bonds_df.iloc[i]['aj'] in atoms_df['atom_id']:
            bonds.append([atoms.Bonds_df.iloc[i]['typ'],\
                atoms.Bonds_df.iloc[i]['ai'], atoms.Bonds_df.iloc[i]['aj']])

# move to mins to orgin
for i in ['x', 'y', 'z']:
    atoms_df[i] = atoms_df[i] - atoms_df[i].min()
# atoms_df['x'] = atoms_df['x'] - atoms_df['x'].min()
print(atoms.Atoms_df.iloc[1]['mol'])
# with concurrent.futures.ThreadPoolExecutor() as executer:
# print(atoms_df['x'].min())
    # executer.map(slice_df, range(ole.NBonds))
end = time.time()

bonds_df = pd.DataFrame(bonds, columns=['typ', 'ai', 'aj'])
bonds_df.index +=1

with open('silica_ole_slic.data', 'w') as f:
    f.write(f"datafile for silica from '{INFILE}'\n")
    f.write(f"\n")
    f.write(f"{len(atoms_df)} atoms\n")
    f.write(f"{max(atoms_df['typ'])} atom types\n")
    f.write(f"{len(bonds_df)} bonds\n")
    f.write(f"{max(bonds_df['typ'])} bond types\n")
    f.write(f"\n")
    f.write(f"{atoms_df['x'].min()} {atoms_df['x'].max()} xlo xhi\n")
    # f.write(f"{ole.Xlim[0]} {ole.Xlim[1]} xlo xhi\n")
    f.write(f"{atoms_df['y'].min()} {atoms_df['y'].max()} ylo yhi\n")
    # f.write(f"{ole.Ylim[0]} {ole.Ylim[1]} ylo yhi\n")
    f.write(f"{atoms_df['z'].min()} {atoms_df['z'].max()} zlo zhi\n")
    # f.write(f"{ole.Zlim[0]} {ole.Zlim[1]} zlo zhi\n")
    f.write(f"\n")
    f.write(f"\n")
    f.write(f"Atoms # full\n")
    f.write(f"\n")
    atoms_df.to_csv(f, index=False, header=None, sep='\t')
    f.write(f"\n")
    f.write(f"Bonds\n")
    f.write(f"\n")
    bonds_df.to_csv(f, index=True, header=None, sep='\t')


print(end - start)
