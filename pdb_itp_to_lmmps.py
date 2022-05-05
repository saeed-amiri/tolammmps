import os, sys, re
import pandas as pd
import numpy as np
import itertools

class DOC:
    """
    reding pdb and itp files from Massimo Delle Piane structure for SiO2 nanoparticles
    """

#GLOBAL VARIABLES and FILE OPERATIONS

ATOM_MASS = dict(HO=1.0080, OB=15.9994, OH=15.9994, OM=15.9994, OMH=15.9994, OD=15.9994, SI=28.0860, SU=28.0860, SD=28.0860)
ATOM_CHARGE = dict(N=0.000, HO=0.400, OH=-0.800, SI=1.600, OB=-0.800, SD=1.500, OD=-1.000, OM=-0.900, OMH=-0.900, SU=1.200)
ATOM_TYPE = dict(HO=1, OB=2, OH=3, OM=4, OMH=5, OD=6, SI=7, SU=8, SD=9)
BOND_TYPE = dict(SIOB=1, SIOH=2, SIOM=3)


def procces_lines(line, lineLen):
        line =line.strip()
        line = line.split(' ')
        line = [item for item in line if item]
        if len(line) != lineLen :
            exit(f'WRONG LINE in line: {line}, EXIT!')
        return line

def drop_digit(obj):
    return re.sub("\d", "", obj)

class PDB:
    """
    reading pdb file and return the coordinates
    """
    def __init__(self, filename) -> None:
        self.filename = filename

    def read_pdb(self):
        self.atom, self.id, self.name, self.label, self.mol, self.x, self.y, self.z = [], [], [], [], [], [], [], [] 
        self.type, self.charge = [], []
        lineCounter = 0
        with open (self.filename, 'r') as f:
            while True:
                line = f.readline()
                lineCounter += 1
                if line.strip().startswith('ATOM'):
                    atom, id, name, label, _, mol, x, y, z, _, _ =  procces_lines(line, 11)
                    name = drop_digit(name)
                    self.atom.append(atom); self.id.append(int(id)); self.name.append(name); self.label.append(label)
                    self.mol.append(mol); self.x.append(float(x)); self.y.append(float(y)); self.z.append(float(z))
                    self.type.append(ATOM_TYPE[name]); self.charge.append(0.0)
                if not line: break
        self.dc = self.make_dict()        
        self.df = self.make_df()
    
    def make_dict(self):
        dc = {'id':self.id, 'mol':self.mol, 'charge':self.charge, 'type':self.type, 'x':self.x, 'y':self.y, 'z':self.z }
        return dc

    def make_df(self):
        return pd.DataFrame.from_dict(self.dc)


class ITP:
    """
    reading itp file to return bonds and angles
    """
    def __init__(self, filename) -> None:
        self.filename = filename

    def read_itp(self):
        atoms, bonds, angles, dihedrals = False, False, False, False
        allAtoms, allBonds, allAngles = [], [], []
        with open(self.filename, 'r') as f:
            while True:
                line = f.readline()
                if line.startswith(';'): continue
                if line.startswith('[ atoms ]'):  atoms, bonds, angles, dihedrals = True, False, False, False
                if line.startswith('[ bonds ]'):  atoms, bonds, angles, dihedrals = False, True, False, False
                if line.startswith('[ angles ]'): atoms, bonds, angles, dihedrals = False, False, True, False
                if line.startswith('[ dihedrals ]'): atoms, bonds, angles, dihedrals = False, False, False, True
                if line.strip() and atoms and not line.startswith('[ '):
                    allAtoms.append(procces_lines(line, 8))
                if line.strip() and bonds and not line.startswith('[ '):
                    allBonds.append(procces_lines(line, 6))
                if line.strip() and angles and not line.startswith('[ '):
                    allAngles.append(procces_lines(line, 8))
                if not line : break
                
        self.itpAtomsDf = self.df_atoms(allAtoms)
        self.itpBondsDf = self.df_bonds(allBonds)
        print(self.itpBondsDf)
        self.itpAnglesDf = self.df_angles(allAngles)
        # self.NmBonds, self.TypBonds, self.SetBonds = self.get_bond_types(self.itpBondsDf['ai_name'].tolist(), self.itpBondsDf['aj_name'].tolist())
        # self.NmAngles, self.TypAngles, self.SetAngles = self.get_angles_types(self.df_angles(allAngles)['ai_name'].tolist(), self.df_angles(allAngles)['aj_name'].tolist(), self.df_angles(allAngles)['ak_name'].tolist())
        # print(self.itpAnglesDf)

    def df_atoms(self, allAtoms) -> pd.DataFrame:
        # making df of all atoms list:
        # [ atoms ]
        # nr      type  resnr resid  atom  cgnr   charge     mass
        nr, type, resnr, resid, atom, cgnr, charge, mass = [], [], [], [], [], [], [], [] 
        for item in allAtoms:
            i_nr,i_type,i_resnr,i_resid,i_atom,i_cgnr,i_charge,i_mass = item
            nr.append(i_nr); type.append(i_type); resnr.append(i_resnr); resid.append(i_resid)
            atom.append(i_atom); cgnr.append(i_cgnr); charge.append(i_charge); mass.append(i_mass)
        # making a dictionary form all the lists
        dic = {'nr':nr, 'type':type, 'resnr':resnr, 'resid':resid, 'atom':atom, 'cgnr':cgnr, 'charge':charge, 'mass':mass}

        del nr, type, resnr, resid, atom, cgnr, charge, mass
        return pd.DataFrame.from_dict(dic)
    
    def df_bonds(self, allBonds) -> pd.DataFrame:
        # making df of all bonds list:
        # [ bonds ]
        # ai        aj        fu ; ai_name aj_name
        id, typ, ai, aj, fu, ai_name, aj_name, bond = [], [], [], [], [], [], [], []
        for counter, item in enumerate(allBonds):
            i_ai, i_aj, i_fu, _, i_ai_name, i_aj_name = item
            id.append(counter+1)
            typ.append(None)
            ai.append(i_ai); aj.append(i_aj); fu.append('#') 
            i_ai_name = drop_digit(i_ai_name)
            i_aj_name = drop_digit(i_aj_name)
            bond.append(f'{i_ai_name}_{i_aj_name}')
            # droping the digits from names
            ai_name.append(i_ai_name); aj_name.append(i_aj_name) 
        # making a dictionary form all the lists
        self.get_bond_types(ai, aj)
        # print(bond)
        dic = {'id': id, 'type':typ, 'ai':ai, 'aj':aj, 'fu':fu, 'bond_name':bond}
        del ai, aj, fu, ai_name, aj_name
        return pd.DataFrame.from_dict(dic)

    def df_angles(self, allAngles) -> pd.DataFrame:
        # making datafram from all angles list:
        # [ angles ]
        # ai        aj        ak  fu ; ai_name aj_name ak_name 
        ai, aj, ak, fu, ai_name, aj_name, ak_name = [], [], [], [], [], [], []
        for item in allAngles:
            i_ai, i_aj, i_ak, i_fu, _, i_ai_name, i_aj_name, i_ak_name = item
            ai.append(i_ai); aj.append(i_aj); ak.append(i_ak); fu.append(i_fu)
            # droping the digits from names
            ai_name.append(drop_digit(i_ai_name)); aj_name.append(drop_digit(i_aj_name)); ak_name.append(drop_digit(i_ak_name))
        # making a dictionary form all the lists
        dic = {'ai':ai, 'aj':aj, 'ak':ak, 'fu':fu, 'ai_name':ai_name, 'aj_name':aj_name, 'ak_name':ak_name}
        del ai, aj, ak, fu, ai_name, aj_name, ak_name
        return pd.DataFrame.from_dict(dic)

    def get_bond_types(self,ai, aj):
        # return number bonds, number of type of bonds, and set of bonds
        bond_names = []
        for i, j in zip(ai, aj):
            bond_names.append(f'{i}_{j}')
        set_of_bonds = set(bond_names)
        number_of_bonds = len(bond_names)
        type_of_bonds = len(set_of_bonds)
        return number_of_bonds, type_of_bonds, set_of_bonds
    
    def get_angles_types(self, ai, aj, ak):
        # return number angles, number of type of angles, and set of angles
        angles_names = []
        for i, j, k in zip(ai, aj, ak):
            angles_names.append(f'{i}_{j}_{k}')
        set_of_angles = set(angles_names)
        number_of_angles = len(angles_names)
        type_of_angles = len(set_of_angles)
        return set_of_angles, number_of_angles, type_of_angles


if __name__ == "__main__":
    pdb = PDB('silica_1nm.pdb')
    pdb.read_pdb()
    itp = ITP('silica_1nm.itp')
    itp.read_itp()