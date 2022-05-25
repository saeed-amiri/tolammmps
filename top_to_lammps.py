import re
import numpy as np

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
    """%FLAG POINTERS                                                                  
    %FORMAT(10I8)                                                                   
    6702       7    5350       0       0       0       0       0       0       0
    8468    3062       0       0       0       3       0       0      10       1
       0       0       0       0       0       0       0       1       8       0
       0

    The cards are space-parsed in some cases, i.e. names of atoms, mols, residues
    """
    def __init__(self) -> None:
        print(f"Reading '{TOPFILE}' ... \n")
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
    """

    def __init__(self) -> None:
        top = READTOP()
        top.get_data()
        self.top = top.FLAG
        del top

    def set_attributes(self) -> None:
        self.get_pointers()
        self.get_atom_name()
        self.get_charges()

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
        # since the format is 10I8 change them in to integer 
        self.top['POINTERS']['data'] = [int(item) for item in self.top['POINTERS']['data']]
        # setting the attributes:
        # set them to None
        length = len(self.top['POINTERS']['data'])
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
        length = len(self.top['CHARGE']['data'])
        kele = 18.2223
        if length != self.NATOM: exit(f"NATOM != N of CHARGE: {length}")
        charges = self.top['CHARGE']['data']
        # convert to float
        charges = [float(q) for q in charges]
        # convert to [e] unit
        charges = [q/kele for q in charges]
        # correct the precision
        charges = [np.round(q, decimals=10) for q in charges]
        self.CHARGE = charges
        del charges
        




if __name__== "__main__":
    TOPFILE = "test3.top"
    top = GETTOP()    
    top.set_attributes()
