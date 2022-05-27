import re
import typing
import numpy as np
import pandas as pd
from pprint import pprint

class DOC:
    """"
    Reading the AMBER data file for SiO2 slab and converting to LAMMPS data file
    """

class TOP:
    """
    reading top file for getting the FLAGS and FORMATS and return a dictionary out of them
    it has 40 FLAG (cards) with different formats in FORTAN style written in the next line
    each flag:
    %FLAG POINTERS                                                                  
    %FORMAT(10I8)                                                                   
    """
    def __init__(self) -> None:
        self.FLAG_list, self.FORMAT_list = [], []

    def read_file(self):
        """
        read line by line and subtract the data from them
        """
        with open (TOPFILE, 'r') as f:
            while True:
                line = f.readline()
                if line.startswith('%'):
                    line = line.split('%')[1]
                    if line.startswith('version'): self.get_version(line)
                    elif line.startswith('FLAG'): self.get_flag(line)
                    elif line.startswith('FORMAT'): self.get_format(line)
                if not line: break
        free_dict = [dict() for _ in self.FORMAT_list]
        self.FLAG = dict(zip(self.FLAG_list, free_dict))
        for i, key in enumerate(self.FLAG): self.FLAG[key]['format'] = self.FORMAT_list[i]
        del free_dict
        del self.FLAG_list
        del self.FORMAT_list

    def get_version(self, line) -> str:
        """geting the version of the AMBER, written in top of the file"""
        self.version = line.strip()
        del line

    def get_flag(self, line) -> list:
        """getting the FLAGEs and make list of it"""
        flag = line.split('FLAG')[1].strip()
        self.FLAG_list.append(flag)
        del line

    def get_format(self, line) -> list:
        """getting the FORMATs and make list of it"""
        format = line.split('FORMAT')[1].strip()
        format = re.findall('\((.*?)\)',format)[0]
        self.FORMAT_list.append(format)
        del line

class READTOP(TOP):
    """
    %FLAG POINTERS                                                                  
    %FORMAT(10I8)                                                                   
    6702       7    5350       0       0       0       0       0       0       0
    8468    3062       0       0       0       3       0       0      10       1
       0       0       0       0       0       0       0       1       8       0
       0

    The cards are space-parsed in some cases, i.e. names of atoms, mols, residues
    """
    def __init__(self) -> None:
        print(f"Reading '{TOPFILE}' ...")
        super().__init__()
        self.read_file()

    def get_data(self) -> None:
        self.mk_modify()
        self.read_card()
        self.crct_card()
        self.mk_format()

    def mk_modify(self) -> None:
        """modify the FLAG dict"""
        self.set_flage()
        self.set_data_list()

    def set_flage(self) -> list:
        """set True and False flag to cards' name"""
        for key in self.FLAG: self.FLAG[key]['flag'] = False

    def set_data_list(self) -> list:
        """giving data atribute to the dict"""
        for key in self.FLAG.keys(): self.FLAG[key]['data']=[]

    def read_card(self) -> list:
        """reading data between two flags"""
        with open (TOPFILE, 'r') as f:
            while True:
                line = f.readline()
                # setting flag = True for the flag we hitting
                if line.startswith("%"):
                    line = line.split("%")[1]
                    if line:
                        line = line.strip()
                        if line.startswith("FLAG"):
                            for key in self.FLAG.keys(): self.FLAG[key]['flag'] = False
                            flag = line.split('FLAG')[1].strip()
                            if flag in self.FLAG.keys(): self.FLAG[flag]['flag'] = True
                        elif line.startswith("FORMAT"): pass
                else:
                    # reading data for each card
                    for key in self.FLAG.keys():
                        # append the line to data file with trimming the trailing newline
                        if self.FLAG[key]['flag']:self.FLAG[key]['data'].append(line.rstrip("\n"))
                if not line: break
    
    def crct_card(self) -> None:
        """correcting the data based on format and removing the blanks"""
        for key in self.FLAG.keys():
            self.FLAG[key]['data'] = ''.join(self.FLAG[key]['data'])
            if self.FLAG[key]['format'] == '20a4': self.do_string(key, 4)
            elif self.FLAG[key]['format'] == '10I8': self.do_string(key, 8)
            elif self.FLAG[key]['format'] == '5E16.8': self.do_string(key, 16)
            elif self.FLAG[key]['format'] == '1a80' :self.do_string(key, 80)
            else: exit(f"\n\tUNDEFINED format in {self.FLAG[key]['format']}\n")

    def do_string(self, key, split) -> list:
        """ fixing the data lists with FORTRAN format: 20a4, 10I8, 5E16.8, 1a80 """
        data_list = self.FLAG[key]['data']
        data_list = [ data_list[i:i+split].strip() for i in range(0, len(data_list)-split+1, split) ]
        self.FLAG[key]['data'] = data_list
        del data_list
    
    def mk_format(self) -> None:
        """correcting the data format """
        for key in self.FLAG.keys():
            if self.FLAG[key]['format'] == '10I8': self.do_integer(key)
            elif self.FLAG[key]['format'] == '5E16.8': self.do_float(key)
            else: pass

    def do_integer(self, key) -> list:
        # fixing the integer format of 10I8
        self.FLAG[key]['data'] = [int(x) for x in self.FLAG[key]['data']]
    
    def do_float(self, key) -> list:
        # fixing the float format of 5E16.8
        self.FLAG[key]['data'] = [float(x) for x in self.FLAG[key]['data']]

    


class GETTOP:
    """
    break down the data information read from TOP file
    This class only extract data from TOPFILE which was read by TOP class.
    Extracting data (anlaysing!), atom-type, ... to be used for making LAMMPS input file will be 
    done by DOTOP class
    """

    def __init__(self, top) -> None:
        self.top = top
        del top
        print("  Setting the attributs ...")

    def get_top(self) -> None:
        # set general attributes which has either NATOM or NRES data-length
        self.set_attributes()
        # get genral data about types: masses, charges, ...
        self.get_types()


    def set_attributes(self) -> None:
        self.get_pointers()
        self.get_atom_name()
        self.get_charges()
        self.get_masses()
        self.get_atom_type()
        self.get_residue_label()
        self.get_residue_pointer()
        self.get_LJ_coeff()

    def get_pointers(self) -> int:
        """
        This section contains the information about how many parameters are present
        in all of the sections. There are 31 or 32 integer pointers (NCOPY might not
        be present). The format and names of all of the pointers are listed below,
        followed by a description of each pointer.
        %FLAG POINTERS
        %FORMAT(10I8)
        NATOM NTYPES NBONH MBONA NTHETH MTHETA NPHIH MPHIA NHPARM NPARM
        NNB NRES NBONA NTHETA NPHIA NUMBND NUMANG NPTRA NATYP NPHB
        IFPERT NBPER NGPER NDPER MBPER MGPER MDPER IFBOX NMXRS IFCAP
        NUMEXTRA NCOPY
        i.e.:
        NATOM   :    Number of atoms
        NTYPES  :    Number of distinct Lennard-Jones atom types
        NBONH   :    Number of bonds containing Hydrogen
        MBONA   :    Number of bonds not containing Hydrogen
        NTHETH  :    Number of angles containing Hydrogen
        MTHETA  :    Number of angles not containing Hydrogen
        NPHIH   :    Number of torsions containing Hydrogen
        MPHIA   :    Number of torsions not containing Hydrogen
        NHPARM  :    Not currently used for anything
        NPARM   :    Used to determine if this is a LES-compatible prmtop
        NNB     :    Number of excluded atoms (length of total exclusion list)
        NRES    :    Number of residues
        NBONA   :    MBONA + number of constraint bonds
        NTHETA  :    MTHETA + number of constraint angles
        NPHIA   :    MPHIA + number of constraint torsions
        NUMBND  :    Number of unique bond types
        NUMANG  :    Number of unique angle types
        NPTRA   :    Number of unique torsion types
        NATYP   :    Number of SOLTY terms. Currently unused.
        NPHB    :    Number of distinct 10-12 hydrogen bond pair types 2 IFPERT Set to 1 if topology contains residue perturbation information.
        NBPER   :    Number of perturbed bonds
        NGPER   :    Number of perturbed angles
        NDPER   :    Number of perturbed torsions
        MBPER   :    Number of bonds in which both atoms are being perturbed
        MGPER   :    Number of angles in which all 3 atoms are being perturbed
        MDPER   :    Number of torsions in which all 4 atoms are being perturbed 3 IFBOX Flag indicating whether a periodic box is present. Values can be 0 (no box), 1 (orthorhombic box) or 2 (truncated octahedro
        NMXRS   :    Number of atoms in the largest residue IFCAP Set to 1 if a solvent CAP is being used
        NUMEXTRA:    Number of extra points in the topology file
        NCOPY   :    Number of PIMD slices or number of beads
        https://ambermd.org/prmtop.pdf
        """
        length = len(self.top['POINTERS']['data'])
        # setting the attributes:
        # set them to None
        nones = lambda n: [None for _ in range(n)]
        self.NATOM, self.NTYPES, self.NBONH, self.MBONA, self.NTHETH, self.MTHETA,\
        self.NPHIH, self.MPHIA, self.NHPARM, self.NPARM, self.NNB, self.NRES, self.NBONA,\
        self.NTHETA, self.NPHIA, self.NUMBND, self.NUMANG, self.NPTRA, self.NATYP, self.NPHB,\
        self.IFPERT, self.NBPER, self.NGPER, self.NDPER, self.MBPER, self.MGPER, self.MDPER,\
        self.IFBOX, self.NMXRS, self.IFCAP, self.NUMEXTRA = nones(31)
        # NCOPY may not be present!
        self.NCOPY = nones(1)
        # setting all the data
        self.NATOM, self.NTYPES, self.NBONH, self.MBONA, self.NTHETH, self.MTHETA,\
        self.NPHIH, self.MPHIA, self.NHPARM, self.NPARM, self.NNB, self.NRES, self.NBONA,\
        self.NTHETA, self.NPHIA, self.NUMBND, self.NUMANG, self.NPTRA, self.NATYP, self.NPHB,\
        self.IFPERT, self.NBPER, self.NGPER, self.NDPER, self.MBPER, self.MGPER, self.MDPER,\
        self.IFBOX, self.NMXRS, self.IFCAP, self.NUMEXTRA = self.top['POINTERS']['data'][0:31]
        if length > 31: self.NCOPY = self.top['POINTERS']['data'][31]

    def get_atom_name(self) -> list:
        """
        make a list of the atoms' name
        his section contains the atom name for every atom in the prmtop.
        %FORMAT(20a4) There are NATOM 4-character strings in this section.
        """
        length = len(self.top['ATOM_NAME']['data'])
        if length != self.NATOM: exit(f"NATOM != N of ATOM_NAME: {length}")
        self.ATOM_NAME = self.top['ATOM_NAME']['data']
        
    def get_charges(self) -> list:
        """
        This section contains the charge for every atom in the prmtop. Charges
        are multiplied by 18.2223 (sqre(k_{ele}) where k_{ele} is the electrostatic constant in
        kcal  ̊A mol^{-1} q^{-2}, where q is the charge of an electron).
        %FORMAT(5E16.8)
        There are NATOM floating-point numbers in this section.
        """
        kele = 18.2223
        charges = self.top['CHARGE']['data']
        length = len(charges)
        if length != self.NATOM: exit(f"NATOM != N of CHARGE: {length}")
        # convert to [e] unit
        charges = [q/kele for q in charges]
        # correct the precision
        charges = [np.round(q, decimals=10) for q in charges]
        self.top['CHARGE']['data'] = charges
        del charges
        
    def get_masses(self) -> list:
        """
        This section contains the atomic mass of every atom in g mol^{-1}.
        %FORMAT(5E16.8)
        There are NATOM floating-point numbers in this section.
        """
        masses = self.top['MASS']['data']
        length = len(masses)
        if length != self.NATOM: exit(f"NATOM != N of MASS: {length}")
        self.top['MASS']['data'] = masses
        del masses
    
    def get_atom_type(self) -> list:
        """
        This section contains the Lennard-Jones atom type index. The Lennard-
        Jones potential contains parameters for every pair of atoms in the system.
        %FORMAT(10I8)
        There are NATOM integers in this section.
        """
        atom_type = self.top['ATOM_TYPE_INDEX']['data']
        length = len(atom_type)
        if length != self.NATOM: exit(f"NATOM != N of ATOM_TYPE_INDEX: {length}")
        self.top['ATOM_TYPE_INDEX']['data'] = atom_type
        del atom_type
    
    def get_residue_label(self) -> list:
        """
        This section contains the residue name for every residue in the prmtop.
        Residue names are limited to 4 letters, and might not be whitespace-delimited
        if any residues have 4-letter names.
        %FORMAT(20a4)
        There are NRES 4-character strings in this section
        """
        length = len(self.top['RESIDUE_LABEL']['data'])
        if length != self.NRES: exit(f"NRES != N of RESIDUE_LABEL: {length}")

    def get_residue_pointer(self) -> list:
        """
        This section lists the first atom in each residue.
        %FORMAT(10i8)
        There are NRES integers in this section
        """
        residue = self.top['RESIDUE_POINTER']['data']
        length = len(residue)
        if length != self.NRES: exit(f"NRES != N of RESIDUE_POINTER: {length}")
        self.top['RESIDUE_POINTER']['data'] = residue
        del residue
        
    def get_LJ_coeff(self) -> list:
        """
        This section contains the LJ A and B-coefficients (ai,j, bi,j in Eq. LENNARD JONES) for all pairs of
        distinct LJ types (see sections ATOM TYPE INDEX and NONBONDED PARM INDEX
        above).
        E_{ij} = [a_{ij}/r^{12}] - [b_{ij}/r^{6}]
        %FORMAT(5E16.8)
        There are [NTYPES * (NTYPES + 1)] /2 floating point numbers in this section.
        !!!!!!!!!!!
        The standard LJ equation in LAMMPS is:
        E_{ij} = 4 * epsilon [ (sigma/r)^{12} - (sigma/r)^{6} ]
        """
        acoeffs = self.top['LENNARD_JONES_ACOEF']['data']
        bcoeffs = self.top['LENNARD_JONES_BCOEF']['data']
        alength = len(acoeffs)
        blength = len(acoeffs)
        len_lj = (self.NTYPES * (self.NTYPES + 1))/2
        if len_lj != alength:exit("\n\tWRONG N of LJ A coeffs\n")
        if len_lj != blength:exit("\n\tWRONG N of LJ B coeffs\n")
        self.LJA = acoeffs
        self.LJB = bcoeffs
        self.top['LENNARD_JONES_ACOEF']['data'] = acoeffs
        self.top['LENNARD_JONES_BCOEF']['data'] = bcoeffs
        self.print_info()
    
    def get_types(self) -> list:
        """
        making a DataFrame from the atoms' name, charges, residues
        to extract the information about the number of types, with their symbols and properties
        """
        data_dict = dict()
        for key in self.top.keys():
            if len(self.top[key]['data'])==self.NATOM:
                data_dict[key]=self.top[key]['data']
        # make a dataframe to extract infos
        self.df = self.mk_df(data_dict)
        del data_dict

    def mk_df(self, data_dict) -> pd.DataFrame:
        df = pd.DataFrame.from_dict(data_dict)
        df = df.drop(['EXCLUDED_ATOMS_LIST'], axis=1)
        df.to_csv('df', sep='\t', index=False)
        df = df.groupby('ATOM_NAME', as_index=False).mean()
        df = df.astype({"ATOM_TYPE_INDEX":int})
        return df


    def print_info(self) -> typing.TextIO:
        print(f"\tseeing {self.NATOM}\t atoms")
        print(f"\tseeing {self.NRES}\t residues")
        print(f"\tseeing {len(set(self.ATOM_NAME))}\t atom names")
        print(f"\tseeing {self.NTYPES}\t atom types")
        print(f"\n")


class PDB:
    """
    reading the PDB file and returning the coordinates
    The PDB file is written in a normal format it is NOT in a standard PDB structure
    """
    def __init__(self) -> None:
        print(f"Reading '{PDBFILE}' ...")

    def read_pdb(self) -> list:
        id, name, residue, chain, x, y, z = [], [], [], [], [], [], []
        lineCounter = 0
        with open (PDBFILE, 'r') as f:
            while True:
                line = f.readline()
                lineCounter += 1
                if line.strip().startswith('ATOM'):
                    line = line.strip().split(' ')
                    line = [item for item in line if item]
                    # ATOM      9 Si   Si      9      -6.272  -1.062  -5.117  1.00  0.00
                    _, i_id, i_name, i_residue, i_chain, i_x, i_y, i_z, _, _ = line
                    id.append(i_id); name.append(i_name); 
                    residue.append(i_residue); chain.append(i_chain)
                    x.append(i_x); y.append(i_y); z.append(i_z) 
                if not line: break

        # set the formats
        id = [int(i) for i in id]
        chain = [int(i) for i in chain]
        x = [float(pos) for pos in x]
        y = [float(pos) for pos in y]
        z = [float(pos) for pos in z]
        # make column for comments
        sharp = ['#' for _ in range(len(x))]
        # make column for charges
        q = [float(0) for _ in range(len(x))]
        # make columns for the flags
        nx = [int(0) for _ in range(len(x))]
        ny = [int(0) for _ in range(len(x))]
        nz = [int(0) for _ in range(len(x))]
        # making a DataFrame in LAMMPS format
        self.mk_lmp_df(id, chain, name, q, x, y, z, nx, ny, nz, sharp)
        del id, chain, name, q, x, y, z, nx, ny, nz, sharp
    
    def mk_lmp_df(self, id, chain, name, q, x, y, z, nx, ny, nz, sharp) -> dict:
        """
        making a DataFrame in LAMMPS 'full' atom style:

        id molecule-tag atom-type q x y z nx ny nz

        Since the atom type is not defined yet, the atom name will be used, and later it will be replaced with the atom type from TOPFILE. 
        Also, all q=0.0 will be set later the charges information will write into a different file.
        'chain' here is the same as molecule-tag in LAMMPS
        """
        data_dict = {'id':id, 'chain':chain, 'name':name, 'q':q, 'x':x, 'y':y, 'z':z, 'nx':nx, 'ny':ny, 'nz':nz, 'sharp':sharp, 'symbol':name}
        self.lmp_df = pd.DataFrame(data_dict)
        self.NATOM = len(id)
        self.NRES = np.max(chain)
        self.ATOM_NAMES= list(set(name))
        self.NNAMES = len(self.ATOM_NAMES)
        del id, chain, name, q, x, y, z, nx, ny, nz, sharp, data_dict
        self.print_info()
    
    def print_info(self) -> typing.TextIO:
        print(f"\t seeing {self.NATOM}\t atoms")
        print(f"\t seeing {self.NRES}\t reseidues (molecules)")
        print(f"\t seeing {self.NNAMES}\t atom names")
        print(f"\n")


class LMPDATA:
    """
    Update DataFrame for LAMMPS input:
    data file in full Atom style:
    id mol type q x y z nx ny nz
    """
    def __init__(self, pdb, top, bonds) -> None:
        self.lmp_df = pdb.lmp_df
        self.top = top
        self.bond_df = bonds
        del pdb, top, bonds

    def mk_lmp(self) -> None:
        self.types, self.charges, self.masses = self.get_q_name_mass()
        # there are Na+ in files, remove them and update all the related attributs
        self.data = self.update_df()
        self.ppo_key('Na+')
        # get box
        self.get_box()
        # get the Number of bonds, index from one
        self.do_bonds()
        # write LAMMPS data file (DATAFILE)
        self.write_data()

    def ppo_key(self, popkey) -> None:
        # checking for the Na+
        try:
            # removing Na+
            self.lmp_df = self.lmp_df[self.lmp_df.symbol != popkey]
            self.NATOM = len(self.lmp_df.id)
            self.NTYPES = self.top.NTYPES - 1
            self.NRES = self.lmp_df.chain.max()
            # remove from masses, charges, and types
            self.types.pop(popkey), self.charges.pop(popkey), self.masses.pop(popkey)
        except: 
            print(f"\tthere were no {popkey} to drop!!\n")
    
    def to_orgin(self) -> None:
        # put mins to origin (0,0,0)
        self.lmp_df.x = self.move_to_zero(self.lmp_df.x)
        self.lmp_df.y = self.move_to_zero(self.lmp_df.y)
        self.lmp_df.z = self.move_to_zero(self.lmp_df.z)


    def move_to_zero(self, data) -> list: return data-np.min(data)

    def update_df(self) -> pd.DataFrame:
        """
        Update DataFram from PDB (self.lmp_df) by substituting the type and charge with date
        from DataFrame from TOP (self.top.df)
        """
        # replace the columns in lmp_df (pdb)
        self.lmp_df['name'], self.lmp_df['q'] = self.set_q_name_mass()

    def get_q_name_mass(self) -> dict:
        types = dict()
        charges = dict()
        masses = dict()
        # make a dict from types and names and charges
        for i, name in enumerate(self.top.df['ATOM_NAME']):
            types[name]=self.top.df.iloc[i]['ATOM_TYPE_INDEX']
            charges[name] = self.top.df.iloc[i]['CHARGE']
            masses[name] = self.top.df.iloc[i]['MASS']
        return types, charges, masses


    def set_q_name_mass(self) -> list:
        # make a list with self.types then replace whole column at once
        typ_lst = []
        q_lst = []
        for name in self.lmp_df['name']:
            typ_lst.append(self.types[name])
            q_lst.append(self.charges[name])
        return typ_lst, q_lst
    
    def get_box(self) -> None:
        # finding box limits
        self.xlo = self.lmp_df['x'].min() - OFFSET; self.xhi = self.lmp_df['x'].max() + OFFSET
        self.ylo = self.lmp_df['y'].min() - OFFSET; self.yhi = self.lmp_df['y'].max() + OFFSET
        self.zlo = self.lmp_df['z'].min() - OFFSET; self.zhi = self.lmp_df['z'].max() + OFFSET

    def do_bonds(self) -> None:
        # set attributes for bonnds
        self.NBONDS = len(self.bond_df)
        self.NBTYPES = self.bond_df['type'].max()
        self.bond_df.index += 1

    def write_data(self) -> typing.TextIO:
        with open(DATAFILE, 'w') as f:
            f.write(f"# dtat from: {TOPFILE} and {PDBFILE}\n")
            f.write(f"\n")
            f.write(f"{self.NATOM} atoms\n")
            f.write(f"{self.NTYPES} atom types\n")
            f.write(f"{self.NBONDS} bonds\n")
            f.write(f"{self.NBTYPES} bond types\n")
            f.write(f"\n")
            f.write(f"{self.xlo} {self.xhi} xlo xhi\n")
            f.write(f"{self.ylo} {self.yhi} ylo yhi\n")
            f.write(f"{self.zlo} {self.zhi} zlo zhi\n")
            f.write(f"\n")
            f.write(f"Atoms # full\n")
            f.write(f"\n")
            self.lmp_df.to_csv(f, sep='\t', index=False, header=None, float_format='%g')
            f.write(f"\n")
            f.write(f"Bonds \n")
            f.write(f"\n")
            self.bond_df.to_csv(f, sep='\t', index=True, header=None, columns=['type', 'ai', 'aj'],float_format='%g')
            self.print_info()

    def print_info(self) -> typing.TextIO:
        print(f"Writing '{DATAFILE}' (LAMMPS data file) ... \n")
        print(f"\t writting {self.NATOM}\t atoms")
        print(f"\t writting {self.NRES}\t reseidues (molecules)")
        print(f"\t writting {self.NTYPES}\t atom types")
        print(f"\n")

class LMPBOND:
    """
    Extract the bonds in the system
    related cards to work with:
     * integer:
        - NTYPES Number of distinct Lennard-Jones atom types
        - NBONH Number of bonds containing Hydrogen
        - MBONA Number of bonds not containing Hydrogen
        - NBONA MBONA + number of constraint bond
     * in cards:
        - BOND_FORCE_CONSTANT
        - BOND_EQUIL_VALUE
        - BONDS_INC_HYDROGEN
        - BONDS_WITHOUT_HYDROGEN
    """

    def __init__(self, top) -> None:
        self.top = top
        del top

    def mk_bonds(self) -> None:
        # print(f"BOND_FORCE_CONSTANT: {self.top.top['BOND_FORCE_CONSTANT']['data']}")
        # print(f"BOND_EQUIL_VALUE: {self.top.top['BOND_EQUIL_VALUE']['data']}")
        # print(f"BONDS_WITHOUT_HYDROGEN: {self.top.top['BONDS_WITHOUT_HYDROGEN']['data']}")
        h_bonds = self.get_h_bond()
        self.bonds = self.mk_hbond_df(h_bonds)

    def get_h_bond(self) -> list:
        """
        BONDS INC HYDROGEN
        This section contains a list of every bond in the system in which at least
        one atom is Hydrogen. Each bond is identified by 3 integers—the two
        atoms involved in the bond and the index into the BOND FORCE CONSTANT
        and BOND EQUIL VALUE. For run-time efficiency, the atom indexes are actu-
        ally indexes into a coordinate array, so the actual atom index A is calculated
        from the coordinate array index N by A = N/3 + 1. (N is the value in the
        topology file)
        %FORMAT(10I8)
        There are 3 * NBONH integers in this section.
        """
        h_bonds = self.top.top['BONDS_INC_HYDROGEN']['data']
        h_bonds = [h_bonds[i:i+3] for i in range(0,len(h_bonds),3)]
        h_bonds = [self.crct_index(lst) for lst in h_bonds]
        return h_bonds
    
    def crct_index(self, lst) -> list:
        for i in range(2):
            lst[i] = self.return_index(lst[i])
        return lst

    def return_index(self, x) -> int: return int((x/3)+1)



    def mk_hbond_df(self, h_bonds) -> pd.DataFrame:
        """return datafram from h_bond list"""
        return pd.DataFrame(h_bonds, columns=['ai', 'aj', 'type'])


class LMPPARAM:
    """
    write out parameters about the system such as 
    mass, coeffs, bonds, ...
    """
    def __init__(self, lmp) -> None:
        self.lmp = lmp
        del lmp

    def mk_types(self) -> None:
        self.get_types()
        self.get_masses()
        self.get_charges()
        self.write_parameters()

    def get_types(self) -> None:
        """
        make clean set of the types, becuase NTYPES != NNAMES
        """
        # convert the the format of the types: int -> str
        self.lmp.types = {k: str(v) for k, v in self.lmp.types.items()}
        # swap keys, and values them tuples the new key
        d  = {tuple(v): k for k, v in self.lmp.types.items()}
        # swap the key and value again
        self.lmp.types = {v: k[0] for k, v in d.items()}
        del d

    def get_masses(self) -> None:
        """ Update masses attribute based on the types """
        self.lmp.masses = {k: v for k, v in self.lmp.masses.items() if k in self.lmp.types.keys()}

    def get_charges(self) -> None:
        """ Update charges (NOT needed, bu incase!) attribute based on the types """
        self.lmp.charges = {k: v for k, v in self.lmp.charges.items() if k in self.lmp.types.keys()}


    def write_parameters(self) -> typing.TextIO:
        with open(PARAMFILE, 'w') as f:
            f.write(f"# Parameters for '{DATAFILE}' from '{TOPFILE}' and '{PDBFILE}'\n")
            f.write(f"\n")
            self.write_mass(f)
            self.write_q(f)
            self.write_group(f)

    def write_mass(self, f) -> typing.TextIO:
        f.write(f"# mass of each type\n")
        for key, value in self.lmp.masses.items():
            f.write(f"mass {self.lmp.types[key]}\t{value:.3f}\t# {self.drop_digit(key)}\n")
        f.write("\n")
    
    def write_q(self, f) -> typing.TextIO:
        f.write(f"# charge of each type\n")
        f.write(f"# charges are already set in data file ({DATAFILE}), here added as comments\n")
        f.write(f"CC\n")
        for key, value in self.lmp.charges.items():
            f.write(f"set type {self.lmp.types[key]}\t{value:.3f}\t# {self.drop_digit(key)}\n")
        f.write(f"CC\n")
        f.write(f"\n")

    def write_group(self, f) -> typing.TextIO:
        # define Hydrogen group
        f.write(f"# defined name of the group \n")
        f.write(f"group Hy_silica type {self.lmp.types['H']}\n")
        f.write(f"\n")

    
    def drop_digit(self, obj) -> str: return re.sub("\d", "", obj)

if __name__== "__main__":
    TOPFILE = "test3.top"
    PDBFILE = "test.pdb"
    DATAFILE = "slab.data"
    PARAMFILE = 'parameters.lmp'
    OFFSET = 0
    data = READTOP()
    data.get_data()
    top = GETTOP(data.FLAG)    
    top.get_top()
    # k = 1
    # for i in range(1, top.NTYPES+1):
        # for j in range(1, top.NTYPES+1):
            # ind = (top.NTYPES * (data.FLAG['ATOM_TYPE_INDEX']['data'][i]-1)) + data.FLAG['ATOM_TYPE_INDEX']['data'][j]
            # print(k, i, j,data.FLAG['NONBONDED_PARM_INDEX']['data'][ind])
            # k+=1
    pdb = PDB()
    pdb.read_pdb()
    lmpbond = LMPBOND(top)
    lmpbond.mk_bonds()
    lmpdata = LMPDATA(pdb, top, lmpbond.bonds)
    lmpdata.mk_lmp()
    lmpparam = LMPPARAM(lmpdata)
    lmpparam.mk_types()
    # pprint(data.FLAG['LENNARD_JONES_BCOEF']['data'])
    # print(len(data.FLAG['ATOM_TYPE_INDEX']['data']))
    
