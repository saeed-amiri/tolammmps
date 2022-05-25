from pprint import pprint
import re
import numpy as np
import pandas as pd
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
                    if line.startswith('FLAG'): self.get_flag(line)
                    if line.startswith('FORMAT'): self.get_format(line)
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

    def get_data(self)->None:
        self.mk_modify()
        self.read_card()
        self.crct_card()

    def mk_modify(self)-> None:
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
        """correcting the data format and removing the blanks"""
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


class GETTOP:
    """
    break down the data information read from TOP file
    This class only extract data from TOPFILE which was read by TOP class.
    Extracting data (anlaysing!), atom-type, ... to be used for making LAMMPS input file will be 
    done by DOTOP class
    """

    def __init__(self) -> None:
        top = READTOP()
        top.get_data()
        self.top = top.FLAG
        del top
        print("  Setting the attributs ...")

    def set_attributes(self) -> None:
        self.get_pointers()
        self.get_atom_name()
        self.get_charges()
        self.get_masses()
        self.get_atom_type()
        self.get_residue_label()
        self.get_residue_pointer()
        self.get_LJ_coeff()

    def get_pointers(self)->int:
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
        # since the format is 10I8 change them in to integer 
        self.top['POINTERS']['data'] = [int(item) for item in self.top['POINTERS']['data']]
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
        # convert to float
        charges = [float(q) for q in charges]
        # convert to [e] unit
        charges = [q/kele for q in charges]
        # correct the precision
        charges = [np.round(q, decimals=10) for q in charges]
        self.CHARGE = charges
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
        masses = [float(m) for m in masses]
        self.MASS = masses
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
        atom_type = [int(atom) for atom in atom_type ]
        self.ATOM_TYPE = atom_type
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
        self.RESIDUE_LABLE = self.top['RESIDUE_LABEL']['data']

    def get_residue_pointer(self):
        """
        This section lists the first atom in each residue.
        %FORMAT(10i8)
        There are NRES integers in this section
        """
        residue = self.top['RESIDUE_POINTER']['data']
        length = len(residue)
        if length != self.NRES: exit(f"NRES != N of RESIDUE_POINTER: {length}")
        residue = [int(item) for item in residue]
        self.RESIDUE_POINTER = residue
        del residue
        
    def get_LJ_coeff(self) -> list:
        """
        This section contains the LJ A and B-coefficients (ai,j, bi,j in Eq. LENNARD JONES) for all pairs of
        distinct LJ types (see sections ATOM TYPE INDEX and NONBONDED PARM INDEX
        above).
        %FORMAT(5E16.8)
        There are [NTYPES * (NTYPES + 1)] /2 floating point numbers in this section.
        """
        acoeffs = self.top['LENNARD_JONES_ACOEF']['data']
        bcoeffs = self.top['LENNARD_JONES_BCOEF']['data']
        alength = len(acoeffs)
        blength = len(acoeffs)
        len_lj = (self.NTYPES * (self.NTYPES + 1))/2
        if len_lj != alength:exit("\n\tWRONG N of LJ A coeffs\n")
        if len_lj != blength:exit("\n\tWRONG N of LJ B coeffs\n")
        acoeffs = [float(a) for a in acoeffs]
        bcoeffs = [float(b) for b in bcoeffs]
        self.LJA = acoeffs
        self.LJB = bcoeffs
        del acoeffs, bcoeffs
        self.print_info()
    
    def print_info(self) -> None:
        print(f"\tseeing {self.NATOM}\t atoms")
        print(f"\tseeing {self.NRES}\t residues")
        print(f"\tseeing {len(set(self.ATOM_NAME))}\t atom names")
        print(f"\tseeing {self.NTYPES}\t atom types")
        print(f"\n")

class DOTOP:
    """
    Extracting data for making LAMMPS inputs
    """

    def __init__(self) -> None:
        top = READTOP()
        top.get_data()
        self.top = top.FLAG
        del top
        print("Extracting data from TOPFILE ...\n")

class PDB:
    """
    reading pdb file and return the coordinates
    The pdb file is written in a normal format it is NOT in a standard PDB structure
    """
    def __init__(self) -> None:
        print(f"Reading '{PDBFILE}' ...")

    def read_pdb(self) -> list:
        id, name, residue, chain, x, y, z = [], [], [], [], [], [], []
        type, charge = [], []
        _sharp, _atom_name = [], []
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
        # put mins to origin (0,0,0)
        x = self.move_to_zero(x)
        y = self.move_to_zero(y)
        z = self.move_to_zero(z)
        # make column for comments
        sharp = ['#' for _ in range(len(x))]
        # make columnm for charges
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
        self.PDBATOMS = pd.DataFrame(data_dict)
        self.NATOM = len(id)
        self.NRES = np.max(chain)
        self.ATOM_NAMES= list(set(name))
        self.NNAMES = len(self.ATOM_NAMES)
        del id, chain, name, q, x, y, z, nx, ny, nz, sharp, data_dict
        self.print_info()
    
    def print_info(self) -> None:
        print(f"\t seeing {self.NATOM}\t atoms")
        print(f"\t seeing {self.NRES}\t reseidues (molecules)")
        print(f"\t seeing {self.NNAMES}\t atom names")
        print(f"\n")

    def move_to_zero(self, data) -> list:
        return data-np.min(data)

if __name__== "__main__":
    TOPFILE = "test3.top"
    PDBFILE = "test.pdb"
    top = GETTOP()    
    top.set_attributes()
    pdb = PDB()
    pdb.read_pdb()
