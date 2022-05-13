import os, sys, re
import pandas as pd
import numpy as np

class DOC:
    """
    A script to correct the bond problem in LAMMPS's write_data command and the boundary conditions.
    These scripts read the data file from write_data, check the distance between particles that have a bond between them;
    if the distance between two particles that have shared a bond is bigger than the HALF OF THE BOX SIZE, move the particle,
    with smaller z close to the other one, i.e., just adding the length of the box to the z component.
    The script is wrote for DATAFILE which have hybrid style, i.e.:
    
    Pair Coeffs # lj/cut/coul/long

    1 lj/cut/coul/long 0.1553 3.166
    2 lj/cut/coul/long 0 0
    
    Bond Coeffs # harmonic

    1 harmonic 600 1
    2 harmonic 268 1.526
    
    Angle Coeffs # harmonic

    1 harmonic 75 109.47
    2 harmonic 58.35 112.7

    Dihedral Coeffs # opls

    1 opls 1.3 -0.05 0.2 0
    2 opls 0 0 0.3 0
    
    usages: {sys.argv[0]} system.data
    """

class FILEERROR:
    """
    there is problem in the header of the DATAFILE,
    maybe long header!\n
    """

class HEADER:
    """
    read haeder data of the data file
    """
    def __init__(self) -> None:
        self.atomsLine = 0
        self.atomsLine = self.check_file()
        print(f'number of header lines: {self.atomsLine}\n')
        self.read_header()


    def check_file(self) -> int:
        FILECHECK = False
        MAXHEADER = 1000
        linecount = 0
        with open(DATAFILE, 'r') as f:
            while True:
                linecount += 1
                line = f.readline()
                if line.startswith('Atoms'):
                    FILECHECK = True
                    atomsLine = linecount
                    break
                if linecount > MAXHEADER:
                    err = FILEERROR()
                    exit(err.__doc__)
                if not line:
                    exit("wrong data file\n")
        return atomsLine

    

    def read_header(self):
        """read header to get all the available info"""
        self.Masses, self.PairCoeff, self.BondCoeff, self.AngleCoeff, self.DihedralCoeff=dict(),dict(),dict(),dict(),dict()
        Masses, PairCoeff, BondCoeff, AngleCoeff, DihedralCoeff, Atoms=False, False, False, False, False, False
        linecount = 0
        with open(DATAFILE, 'r') as f:
            while True:
                linecount += 1
                line = f.readline()
                if line.strip().endswith("atoms"):      self.NATOMS = int(line.strip().split(' ')[0])
                if line.strip().endswith("atom types"): self.NATomTyp = int(line.strip().split(' ')[0])
                if line.strip().endswith("bonds"):      self.NBonds = int(line.strip().split(' ')[0])
                if line.strip().endswith("bond types"): self.NBondTyp = int(line.strip().split(' ')[0])
                if line.strip().endswith("angles"):     self.NAngles = int(line.strip().split(' ')[0])
                if line.strip().endswith("angle types"):self.NAngleTyp = int(line.strip().split(' ')[0])
                if line.strip().endswith("dihedrals"):  self.NDihedrals = int(line.strip().split(' ')[0])
                if line.strip().endswith("dihedral typss"): self.NDihedralsTyp = int(line.strip().split(' ')[0])
                if line.strip().endswith("xhi"): self.Xlim = self.get_axis_lim(line.strip().split('xlo')[0])
                if line.strip().endswith("yhi"): self.Ylim = self.get_axis_lim(line.strip().split('ylo')[0])
                if line.strip().endswith("zhi"): self.Zlim = self.get_axis_lim(line.strip().split('zlo')[0])

                if line.strip().startswith("Masses"):
                    Masses, PairCoeff, BondCoeff, AngleCoeff, DihedralCoeff, Atoms = True, False, False, False, False, False
                if line.strip().startswith("Pair"):
                    Masses, PairCoeff, BondCoeff, AngleCoeff, DihedralCoeff, Atoms = False, True, False, False, False, False
                if line.strip().startswith("Bond Coeffs"):
                    Masses, PairCoeff, BondCoeff, AngleCoeff, DihedralCoeff, Atoms = False, False, True, False, False, False
                if line.strip().startswith("Angle Coeffs"):
                    Masses, PairCoeff, BondCoeff, AngleCoeff, DihedralCoeff, Atoms = False, False, False, True, False, False
                if line.strip().startswith("Dihedral Coeffs"):
                    Masses, PairCoeff, BondCoeff, AngleCoeff, DihedralCoeff, Atoms = False, False, False, False, True, False
                if line.strip().startswith("Atoms"): 
                    Masses, PairCoeff, BondCoeff, AngleCoeff, DihedralCoeff, Atoms = False, False, False, False, False, True
                if line.strip():
                    if Masses: self.get_masses(line.strip(), 'Masses')
                    if PairCoeff: self.get_pair_coeff(line.strip(), 'Pair')
                    if BondCoeff: self.get_bond_coeff(line.strip(), 'Bond')
                    if AngleCoeff: self.get_angle_coeff(line.strip(), 'Angle')
                    if DihedralCoeff: self.get_dihedral_coeff(line.strip(), 'Dihedral')
                if Atoms: break
                if not line: break

    def get_axis_lim(self, lim) -> list:
        lim = lim.split(' ')
        lim = [float(item) for item in lim if item]
        return lim
    
    def get_masses(self, line, check) -> dict:
        if check not in line: 
            typ = line.split(' ')[0]
            mass = float(line.split(' ')[1])
            self.Masses[typ]=mass
        else: pass

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




class BODY():
    """
    read the data for atoms,velocities, bonds, angles, dihedrals
    """

    def __init__(self) -> None:
        pass
    
    def read_body(self):
        self.Atoms, self.Velocities, self.Bonds, self.Angles, self.Dihedrals = dict(), dict(), dict(), dict(), dict()
        Atoms, Velocities, Bonds, Angles, Dihedrals = False, False, False, False, False
        with open(DATAFILE, 'r') as f:
            while True:
                line = f.readline()
                if line.strip().startswith('Atoms'):
                    Atoms, Velocities, Bonds, Angles, Dihedrals = True, False, False, False, False
                if line.strip().startswith('Velocities'):
                    Atoms, Velocities, Bonds, Angles, Dihedrals = False, True, False, False, False
                if line.strip().startswith('Bonds'):
                    Atoms, Velocities, Bonds, Angles, Dihedrals = False, False, True, False, False
                if line.strip().startswith('Angles'):
                    Atoms, Velocities, Bonds, Angles, Dihedrals = False, False, False, True, False
                if line.strip().startswith('Dihedrals'):
                    Atoms, Velocities, Bonds, Angles, Dihedrals = False, False, False, False, True
                if line.strip():
                    if Atoms :self.get_atoms(line.strip())
                    if Velocities: self.get_velocities(line.strip())
                    if Bonds: self.get_bonds(line.strip())
                    if Angles: self.get_angles(line.strip())
                    if Dihedrals: self.get_dihedrals(line.strip())
                if not line: break
    
    def get_atoms(self, line) -> dict:
        if 'Atoms' not in line:
            line = line.split(" ")
            line = [item for item in line if item]
            atom_id = int(line[0])
            i_mol=int(line[1]); i_typ=int(line[2]); i_charg=float(line[3]) 
            i_x=float(line[4]);i_y=float(line[5]); i_z=float(line[6])
            i_nx=str(line[7]); i_ny=str(line[8]); i_nz=str(line[9])
            self.Atoms[atom_id]=dict(mol=i_mol, typ=i_typ, charg=i_charg, x=i_x,y=i_y, z=i_z, nx=i_nx, ny=i_ny, nz=i_nz)
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
            self.Angles[angle_id] = dict(typ = i_typ, ai = i_ai, aj = i_aj, ak = i_ak)

    def get_dihedrals(self, line) -> dict:
        if "Dihedrals" not in line:
            line = line.split(" ")
            line = [int(item) for item in line if item]
            dihedrals_id = line[0]
            i_typ = line[1]; i_ai = line[2]; i_aj = line[3]; i_ak = line[4]; i_ah = line[5]
            self.Dihedrals[dihedrals_id] = dict(typ = i_typ, ai = i_ai, aj = i_aj, ak = i_ak, ah = i_ah)

class UPDATE:
    """
    Update data, checking the bonds
    """
    def __init__(self, Bonds, Atoms, header) -> None:
        self.Bonds = Bonds
        self.Atoms = Atoms
        self.header = header
        self.get_sizes()
        del Bonds, Atoms, header

    def get_sizes(self):
        self.boxx = self.header.Xlim[1]-self.header.Xlim[0]
        self.boxy = self.header.Ylim[1]-self.header.Ylim[0]
        self.boxz = self.header.Zlim[1]-self.header.Zlim[0]
        print(self.boxx, self.boxy, self.boxz)

    def check_bonds(self):
        for bond in range(1,len(self.Bonds)+1):
            ai = self.Bonds[bond]['ai']
            aj = self.Bonds[bond]['aj']
            self.check_coords(ai, aj)

    def check_coords(self, ai, aj):
        axis = ['x', 'y', 'z']
        lims = [self.boxx, self.boxy, self.boxz]
        for ax, lim in zip(axis, lims):
            xi = self.Atoms[ai][ax]
            xj = self.Atoms[aj][ax]
            distance = np.abs(xi - xj)
            if self.Atoms[ai]['mol'] != self.Atoms[aj]['mol']:
                exit(f'DIFFRENT MOLECULE ID{ai}-{aj}')
            if distance > lim*0.25:
                if xi < xj: xi += lim
                else: xj += lim
                self.Atoms[ai][ax] = xi
                self.Atoms[aj][ax] = xj

if __name__ == "__main__":
    # check the input file 
    if len(sys.argv)==1:
        doc = DOC()
        exit(f'\nONE INPUT IS RWUIRED\n{doc.__doc__}')
    DATAFILE = sys.argv[1]
    OUTEX = DATAFILE.split('.')[0]
    header = HEADER()
    body = BODY()
    body.read_body()
    df = pd.DataFrame(body.Atoms).T
    df = df.sort_index(axis=0)
    df.to_csv(f'atoms.{OUTEX}', index=True, header=None, sep=' ')
    update = UPDATE(body.Bonds, body.Atoms, header)
    update.check_bonds()
    df = pd.DataFrame(update.Atoms).T
    df = df.sort_index(axis=0)
    df.to_csv(f'bonds.{OUTEX}', index=True, header=None, sep=' ')

