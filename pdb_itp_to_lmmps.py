from concurrent.futures import thread
import itertools
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


def procces_lines(line, lineLen) -> list:
        line =line.strip()
        line = line.split(' ')
        line = [item for item in line if item]
        if len(line) != lineLen :
            exit(f'WRONG LINE in line: {line}, EXIT!')
        return line
    

def drop_digit(obj) -> str:
    return re.sub("\d", "", obj)

def drop_semicolon(line)->list:
    return re.sub(r'\;.*', "", line)

def return_df_value(df, name, name_i,  string_re) -> float:
    # df: dataframe
    # name: name of the column to look for a match
    # name_i: the string to look for in the 'name' column
    # string_re: return the value of the selected coulmn
    return df.loc[df[name]==name_i][string_re].values[0]

def mix_epsilon(si, sj) -> float: return np.sqrt(si * sj)
def mix_sigma_geometric(si, sj) -> float: return np.sqrt(si * sj)
def mix_sigma_arithmetic(si, sj) -> float: return 0.5*(si + sj)

def move_to_zero(data) -> list:
    return data-np.min(data)

class PDB:
    """
    reading pdb file and return the coordinates
    """
    def __init__(self, filename) -> None:
        self.filename = filename

    def read_pdb(self) -> list:
        self.atom, self.id, self.name, self.label, self.mol, self.x, self.y, self.z = [], [], [], [], [], [], [], [] 
        self.type, self.charge = [], []
        self._sharp, self._atom_name = [], []
        lineCounter = 0
        with open (self.filename, 'r') as f:
            while True:
                line = f.readline()
                lineCounter += 1
                if line.strip().startswith('ATOM'):
                    # atom_serial,atom_name, chain_id, residue_number, residue_name, x, y, z
                    id, name, chain_id, mol, _, x, y, z = self.procces_pdb(line)
                    name = drop_digit(name)
                    self.id.append(int(id)); self.name.append(name); self.label.append(chain_id)
                    self.mol.append(mol); self.x.append(float(x)); self.y.append(float(y)); self.z.append(float(z))
                    self.type.append(ATOM_TYPE[name]); self.charge.append(0.0)
                    self._sharp.append('#'); self._atom_name.append(name)
                if not line: break
        self.x = move_to_zero(self.x)
        self.y = move_to_zero(self.y)
        self.z = move_to_zero(self.z)
        self.dc = self.make_dict()        
        self.df = self.make_df()
        self.NAtoms = len(self.id)
        self.xlo = np.min(self.x); self.xhi = np.max(self.x)
        self.ylo = np.min(self.y); self.yhi = np.max(self.y)
        self.zlo = np.min(self.z); self.zhi = np.max(self.z)

    def procces_pdb(self, line) -> list:
        atom_serial = line[6:11].strip()
        atom_name = line[12:16].strip()
        residue_name = line[17:20].strip()
        chain_id = line[21].strip()
        residue_number = line[22:26].strip()
        x = line[30:38].strip()
        y = line[38:46].strip()
        z = line[46:54].strip()
        return atom_serial,atom_name, chain_id, residue_number, residue_name, x, y, z
    
    def make_dict(self) -> dict:
        dc = {'#id':self.id, 'mol':self.mol, 'type':self.type, 'charge':self.charge, \
             'x':self.x, 'y':self.y, 'z':self.z, 'com':self._sharp, 'name':self._atom_name}
        return dc

    def make_df(self) -> pd.DataFrame:
        return pd.DataFrame.from_dict(self.dc)


class ITP:
    """
    reading itp file to return bonds and angles
    """
    def __init__(self, filename) -> None:
        self.filename = filename

    def read_itp(self) -> pd.DataFrame:
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
        self.itpBondsDf = self.make_df_bonds(allBonds)
        self.itpAnglesDf = self.make_df_angles(allAngles)

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
    
    def make_df_bonds(self, allBonds) -> pd.DataFrame:
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
            ai_name.append(i_ai_name)
            aj_name.append(i_aj_name) 

        # getting the number bonds infos
        self.NmBonds, self.TypBonds, self.SetBonds = self.get_bond_types(bond)
        # updating the type of eche bond 
        for i, (_, n) in enumerate(zip(typ,bond)):
            typ[i] = self.SetBonds[n]
        # making a dictionary form all the lists
        dic = {'id': id, 'type':typ, 'ai':ai, 'aj':aj, 'fu':fu, 'bond_name':bond}
        
        del ai, aj, fu, ai_name, aj_name
        return pd.DataFrame.from_dict(dic)

    def make_df_angles(self, allAngles) -> pd.DataFrame:
        # making datafram from all angles list:
        # [ angles ]
        # ai        aj        ak  fu ; ai_name aj_name ak_name 
        id, typ, ai, aj, ak, fu, ai_name, aj_name, ak_name, angle = [], [], [], [], [], [], [], [], [], []
        for counter, item in enumerate(allAngles):
            i_ai, i_aj, i_ak, i_fu, _, i_ai_name, i_aj_name, i_ak_name = item
            id.append(counter+1)
            typ.append(None)
            ai.append(i_ai); aj.append(i_aj); ak.append(i_ak); fu.append("#")
            # droping the digits from names
            i_ai_name = drop_digit(i_ai_name)
            i_aj_name = drop_digit(i_aj_name)
            i_ak_name = drop_digit(i_ak_name)
            ai_name.append(i_ai_name)
            aj_name.append(i_aj_name)
            ak_name.append(i_ak_name)
            angle.append(f'{i_ai_name}_{i_aj_name}_{i_ak_name}')
        self.NmAngles, self.TypAngles, self.SetAngles = self.get_angles_types(angle)
        # updating the type of eche angel 
        for i, (_, n) in enumerate(zip(typ,angle)):
            typ[i] = self.SetAngles[n]
        # making a dictionary form all the lists
        dic = {'id':id, 'typ':typ, 'ai':ai, 'aj':aj, 'ak':ak, 'fu':fu, 'angle_name':angle}
        
        del ai, aj, ak, fu, ai_name, aj_name, ak_name
        return pd.DataFrame.from_dict(dic)

    def get_bond_types(self,bond_names) -> list:
        # return number bonds, number of type of bonds, and set of bonds
        set_of_bonds = set(bond_names)
        number_of_bonds = len(bond_names)
        type_of_bonds = len(set_of_bonds)
        # convert set of bonds to dictionary
        set_of_bonds = {item:i+1 for i,item in enumerate(set_of_bonds)}
        return number_of_bonds, type_of_bonds, set_of_bonds
    
    def get_angles_types(self, angles_names) -> list:
        # return number angles, number of type of angles, and set of angles
        set_of_angles = set(angles_names)
        number_of_angles = len(angles_names)
        type_of_angles = len(set_of_angles)
        # convert set of angles to dictionary
        set_of_angles = {item:i+1 for i,item in enumerate(set_of_angles)}
        return number_of_angles, type_of_angles, set_of_angles

class CHARMM:
    """
    reading charmm36_silica.itp
    to get interactions parameters
    """

    def __init__(self) -> None:
        pass

    def read_charmm(self) -> pd.DataFrame:
        atomtypes, nonbond_params, bondtypes, pairtypes, angletypes, dihedraltypes = False, False, False, False, False, False
        atomList, nonbondList, bondList, pairList,  angleList, dihedralList = [], [], [], [], [], []
        with open(CHARMMFILE, 'r') as f:
            while True:
                line = f.readline()
                if line.startswith(';'): continue
                if line.startswith('[ atomtypes ]'): 
                    atomtypes, nonbond_params, bondtypes, pairtypes, angletypes, dihedraltypes = True, False, False, False, False, False
                if line.startswith('[ nonbond_params ]'): 
                    atomtypes, nonbond_params, bondtypes, pairtypes, angletypes, dihedraltypes = False, True, False, False, False, False
                if line.startswith('[ bondtypes ]'): 
                    atomtypes, nonbond_params, bondtypes, pairtypes, angletypes, dihedraltypes = False, False, True, False, False, False
                if line.startswith('[ pairtypes ]'): 
                    atomtypes, nonbond_params, bondtypes, pairtypes, angletypes, dihedraltypes = False, False, False, True, False, False
                if line.startswith('[ angletypes ]'): 
                    atomtypes, nonbond_params, bondtypes, pairtypes, angletypes, dihedraltypes = False, False, False, False, True, False
                if line.startswith('[ dihedraltypes ]'): 
                    atomtypes, nonbond_params, bondtypes, pairtypes, angletypes, dihedraltypes = False, False, False, False, False, True
                if line.strip() and atomtypes and not line.startswith('[ '):
                    line = drop_semicolon(line).strip()
                    atomList.append(line)
                if line.strip() and nonbond_params and not line.startswith('[ '):
                    nonbondList.append(line)
                if line.strip() and bondtypes and not line.startswith('[ '):
                    bondList.append(line)
                if line.strip() and pairtypes and not line.startswith('[ '):
                    pairList.append(line)
                if line.strip() and angletypes and not line.startswith('[ '):
                    angleList.append(line)
                if line.strip() and dihedraltypes and not line.startswith('[ '):
                    dihedralList.append(line)
                if not line: break
        self.atomtyps_df = self.read_atomtypes(atomList); del atomList
        self.nonbond_df = self.read_nonbond(nonbondList); del nonbondList
        self.bond_df = self.read_bondtypes(bondList); del bondList
        self.pair_df = self.read_pairtypes(pairList); del pairList
        self.angle_df = self.read_angletypes(angleList); del angleList


    def read_atomtypes(self, atomList) -> pd.DataFrame:
        atom_dict = dict(name=[] ,at_num=[] ,mass=[] ,charge=[] ,ptype=[] ,sigma=[] ,epsilon=[])
        for item in atomList:
            i_name, i_at_num, i_mass, i_charge, i_ptype, i_sigma, i_epsilon = procces_lines(item,7)
            i_name = i_name.strip()
            i_at_num = int(i_at_num.strip())
            i_mass = float(i_mass.strip())
            i_charge = float(i_charge.strip())
            i_ptype = i_ptype.strip()
            i_sigma = float(i_sigma.strip())
            i_epsilon = float(i_epsilon.strip())
            atom_dict['name'].append(i_name); atom_dict['at_num'].append(i_at_num)
            atom_dict['mass'].append(i_mass); atom_dict['charge'].append(i_charge)
            atom_dict['ptype'].append(i_ptype); atom_dict['sigma'].append(i_sigma)
            atom_dict['epsilon'].append(i_epsilon)
        del atomList
        return pd.DataFrame.from_dict(atom_dict)    
    
    def read_nonbond(self, nonbondList) -> pd.DataFrame:
        nonbond_dict = dict(ai=[], aj=[], func=[], sigma=[], epsilon=[])
        for item in nonbondList:
            i_ai, i_aj, i_func, i_sigma, i_epsilon, = procces_lines(item, 5)
            i_ai = i_ai.strip()
            i_aj = i_aj.strip()
            i_func = i_func.strip()
            i_sigma = float(i_sigma.strip())
            i_epsilon = float(i_epsilon.strip())
            nonbond_dict['ai'].append(i_ai); nonbond_dict['aj'].append(i_aj)
            nonbond_dict['func'].append(i_func); nonbond_dict['sigma'].append(i_sigma)
            nonbond_dict['epsilon'].append(i_epsilon)
        del nonbondList
        return pd.DataFrame.from_dict(nonbond_dict)


        
    def read_bondtypes(self, bondList) -> pd.DataFrame:
        bond_dict = dict(bond=[], func=[], b0=[], Kb=[])
        for item in bondList:
            i_ai, i_aj, i_func, i_b0, i_Kb = procces_lines(item, 5)
            i_ai = i_ai.strip()
            i_aj = i_aj.strip()
            i_func = i_func.strip()
            i_b0 = float(i_b0.strip())
            i_Kb = float(i_Kb.strip())
            bond_dict['bond'].append(f'{i_ai}_{i_aj}')
            bond_dict['func'].append(i_func); bond_dict['b0'].append(i_b0)
            bond_dict['Kb'].append(i_Kb)
        del bondList
        return pd.DataFrame.from_dict(bond_dict)

    def read_pairtypes(self, pairList) -> pd.DataFrame:
        pair_dict = dict(pairs=[], func=[], sigma=[], epsilon=[])
        for item in pairList:
            i_ai, i_aj, i_func, i_sigma, i_epsilon = procces_lines(item, 5)
            i_ai = i_ai.strip()
            i_aj = i_aj.strip()
            i_func = i_func.strip()
            i_sigma = float(i_sigma.strip())
            i_epsilon = float(i_epsilon.strip())
            i_pairs = f'{i_ai}_{i_aj}'
            pair_dict['pairs'].append(i_pairs)
            pair_dict['func'].append(i_func); pair_dict['sigma'].append(i_sigma)
            pair_dict['epsilon'].append(i_epsilon)
        del pairList
        return pd.DataFrame.from_dict(pair_dict)

    def read_angletypes(self, angleList) -> pd.DataFrame:
        angle_dict = dict(angle=[], func=[], th0=[], cth=[], S0=[], Kub=[])
        for item in angleList:
            i_ai, i_aj, i_ak, i_func, i_th0, i_cth, i_S0, i_Kub = procces_lines(item, 8)
            i_ai = i_ai.strip()
            i_aj = i_aj.strip()
            i_ak = i_ak.strip()
            i_func = i_func.strip()
            i_th0 = float(i_th0.strip())
            i_cth = float(i_cth.strip())
            i_S0 = float(i_S0.strip())
            i_Kub = float(i_Kub.strip())
            i_angle = f'{i_ai}_{i_aj}_{i_ak}'
            angle_dict['angle'].append(i_angle)
            angle_dict['func'].append(i_func); angle_dict['th0'].append(i_th0); angle_dict['cth'].append(i_cth)
            angle_dict['S0'].append(i_S0); angle_dict['Kub'].append(i_Kub)
        del angleList
        return pd.DataFrame.from_dict(angle_dict)

    def read_dihedraltypes(self):
        pass
                


class WRITE_DATA :
    """
    write out the output file for lammps
    """
    def __init__(self, pdb, itp) -> None:
        self.pdb = pdb        
        self.itp = itp
        del pdb, itp

    def write_file(self) -> None:
        with open(DATAFILE, 'w') as f:
            f.write(f'#\twrite data file with "{sys.argv[0]}" form {ITPFILE} and {PDBFILE}\n')
            f.write(f'#\tmasses, interactions (bond and non-bond) are written in {PARAMFILE}\n\n')
            self.write_numbers_types(f)
            self.write_box(f)
            self.write_masses(f)
            self.write_pair_coef(f)
            self.write_bonds_coef(f)
            self.write_coords(f)
            self.write_bonds(f)
            self.write_angles(f)

    
    def write_numbers_types(self, f) -> None:
        f.write(f'{self.pdb.NAtoms}\tatoms\n')
        f.write(f'{len(ATOM_TYPE)}\tatom types\n')
        f.write(f'{self.itp.NmBonds}\tbonds\n')
        f.write(f'{itp.TypBonds}\tbond types\n')
        f.write(f'{self.itp.NmAngles}\tangles\n')        
        f.write(f'{itp.TypAngles}\tangle types\n')
        f.write(f'\n')

    def write_box(self, f):
        f.write(f'{self.pdb.xlo}\t{self.pdb.xhi}\txlo\txhi\n')
        f.write(f'{self.pdb.ylo}\t{self.pdb.yhi}\tylo\tyhi\n')
        f.write(f'{self.pdb.zlo}\t{self.pdb.zhi}\tzlo\tzhi\n')
        f.write(f'\n')

    def write_masses(self, f) -> None:
        f.write(f'Masses\n\n')
        for i, (m, t) in enumerate(zip(ATOM_MASS, ATOM_TYPE)):
            if m != t : exit(f"WRONG mass and type sets: {t} :: {m}")
            f.write(f'{ATOM_TYPE[m]}\t{ATOM_MASS[t]}\t# {m}\n')
        f.write(f'\n')    

    def write_pair_coef(self, f) -> None:
        pairs = self.make_pairs()
        f.write(f'PairIJ Coeffs\n\n')
        for i, item in enumerate(list(pairs)):
            ai = [name for name, id in ATOM_TYPE.items() if id == item[0]][0]
            aj = [name for name, id in ATOM_TYPE.items() if id == item[1]][0]
            _pairIJ = f'{ai}_{aj}'
            try:
                sigma = return_df_value(charmm.pair_df, 'pairs', _pairIJ, 'sigma' )
                epsilon = return_df_value(charmm.pair_df, 'pairs', _pairIJ, 'epsilon' )
            except:
                sigma = 0.0
                epsilon = 1.0
            if _pairIJ in self.itp.SetBonds: exit('EXIT!! there is a pair with bonding and non-bonding interactions!!')
            f.write(f'{item[0]}\t{item[1]}  {epsilon} {sigma} # {_pairIJ} \n')
        f.write('\n')

    def write_bonds_coef(self, f) -> None:
        f.write('Bonds Coeffs # harmonic\n\n')
        for i, item in enumerate(self.itp.SetBonds):
            Kb = return_df_value(charmm.bond_df, 'bond', item, 'Kb')
            b0 = return_df_value(charmm.bond_df, 'bond', item, 'b0')
            f.write(f'{i+1} {Kb} {b0} # {item}\n')
        f.write('\n')


    def write_coords(self, f) -> None:
        f.write(f'Atoms\t#\tfull\n\n')
        self.pdb.df.to_csv(f, mode='a',header=None, index=False, sep=' ')
        f.write(f'\n')

    def write_bonds(self, f) -> None:
        f.write(f'Bonds\n\n')
        self.itp.itpBondsDf.to_csv(f, mode='a', header=None, index=False, sep=' ')
        f.write(f'\n')

    def write_angles(self, f) -> None:
        f.write(f'Angles\n\n')
        self.itp.itpAnglesDf.to_csv(f, mode='a', header=None, index=False, sep=' ')

    def make_pairs(self) -> None:
        _type_list = [i+1 for i, _ in enumerate(ATOM_TYPE) ]
        return itertools.combinations_with_replacement(_type_list, 2)




if len(sys.argv) < 2:
    err = DOC()
    print(err.__doc__)
    exit(f'USAGE Ex.:\n {sys.argv[0]} silica_1nm\n')

SIO2 = sys.argv[1]
ITPFILE = f'{SIO2}.itp'
PDBFILE = f'{SIO2}.pdb'
DATAFILE = f'{SIO2}.data'
PARAMFILE = f'{SIO2}.data'
CHARMMFILE = 'charmm36_silica.itp'

if __name__ == "__main__":
    # MAKE GLOBAL PARAMETERS FROM CHARMM DATA FILE
    charmm = CHARMM()
    charmm.read_charmm()
    itp = ITP(ITPFILE)
    itp.read_itp()
    pdb = PDB(PDBFILE)
    pdb.read_pdb()
    # param = WRITE_PARAM(pdb, itp)
    # param.write_param()
    out = WRITE_DATA(pdb, itp)
    out.write_file()
