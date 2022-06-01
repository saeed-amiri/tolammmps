import itertools
import re
import typing
import numpy as np
import pandas as pd
from pprint import pprint


class DOC:
    """" Converting data:
    Reading the AMBER data file for SiO2 slab and converting to LAMMPS data
    file.
    The script read data for the Silica slab. in a pdb and top (AMBER) format.
    The script extract data based on the these refrenceses:
        - https://ambermd.org/FileFormats.php#topo.cntrl
        - https://ambermd.org/prmtop.pdf

    and write out the output based for the LAMMPS, based on the:
        - https://ambermd.org/prmtop.pdf

    """


class TOP:
    """reading top file for getting the FLAGS and FORMATS and return a
    dictionary out of them it has 40 FLAG (cards) with different formats
    in FORTAN style written in the next line for each flag:
    %FLAG POINTERS
    %FORMAT(10I8)
    """
    def __init__(self) -> None:
        self.FLAG_list, self.FORMAT_list = [], []

    def read_file(self) -> None:
        """ read line by line and subtract the data from them """
        with open(TOPFILE, 'r') as f:
            while True:
                line = f.readline()
                self.process_line(line)
                if not line:
                    break
        # making a list of dictionaries to save data of each flag
        # we do not know the numbers of flags yet
        empty_dict = [dict() for _ in self.FORMAT_list]
        # set atribeuts for the for all the falgs
        self.FLAG = dict(zip(self.FLAG_list, empty_dict))
        # append the format of each falg to the dictionary if esch one of them
        for i, key in enumerate(self.FLAG):
            self.FLAG[key]['format'] = self.FORMAT_list[i]
        # flush the memory
        del empty_dict
        del self.FLAG_list
        del self.FORMAT_list

    def process_line(self, line) -> None:
        # process line for appending the data to lists and dicts
        if line.startswith('%'):
            line = line.split('%')[1]
            if line.startswith('version'):
                self.get_version(line)
            elif line.startswith('FLAG'):
                self.get_flag(line)
            elif line.startswith('FORMAT'):
                self.get_format(line)
        else: pass

    def get_version(self, line: list) -> str:
        """geting the version of the AMBER, written in top of the file"""
        self.version = line.strip()
        del line

    def get_flag(self, line: list) -> list:
        """getting the FLAGEs and make list of it"""
        flag = line.split('FLAG')[1].strip()
        self.FLAG_list.append(flag)
        del line

    def get_format(self, line: list) -> list:
        """getting the FORMATs and make list of it"""
        format = line.split('FORMAT')[1].strip()
        format = re.findall('\((.*?)\)', format)[0]
        self.FORMAT_list.append(format)
        del line


class READTOP(TOP):
    """Rendering data from TOP class.
    Any FALG needed later on should proceses through this class
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
        # call functions one by one to get data
        # These functions return None; they just set attributes
        self.mk_modify()
        self.read_card()
        self.crct_card()
        self.mk_format()

    def mk_modify(self) -> None:
        """modify the FLAG dict"""
        self.set_flage()
        self.set_data_list()

    def set_flage(self) -> None:
        """set True and False flag to cards' name"""
        for key in self.FLAG:
            self.FLAG[key]['flag'] = False

    def set_data_list(self) -> None:
        """giving data atribute to the dict"""
        for key in self.FLAG.keys():
            self.FLAG[key]['data'] = []

    def read_card(self) -> None:
        """reading data between two flags
        class names:
            FLAG: using for each FLAG header in data file
            falg (variable): real True and False falg to indicate
            which data is reading
            'flag': using for save True and False state 
        """
        with open(TOPFILE, 'r') as f:
            while True:
                line = f.readline()
                # setting flag = True for the flag we hitting
                if line.startswith("%"):
                    line = line.split("%")[1]
                    self.process_flag(line)
                else:
                    # reading data for each card
                    self.grap_data(line)
                if not line:
                    break
    
    def process_flag(self, line: list) -> None:
        if line:
            line = line.strip()
            if line.startswith("FLAG"):
                for key in self.FLAG.keys():
                    self.FLAG[key]['flag'] = False
                flag = line.split('FLAG')[1].strip()
                if flag in self.FLAG.keys():
                    self.FLAG[flag]['flag'] = True
            elif line.startswith("FORMAT"):
                pass
    
    def grap_data(self, line: list) -> None:
        # append data for each FLAG:
        for key in self.FLAG.keys():
            # append the line to data file with trimming the trailing newline
            if self.FLAG[key]['flag']:
                self.FLAG[key]['data'].append(line.rstrip("\n"))

    def crct_card(self) -> None:
        """correcting the data based on format and removing the blanks"""
        for key in self.FLAG.keys():
            self.FLAG[key]['data'] = ''.join(self.FLAG[key]['data'])
            if self.FLAG[key]['format'] == '20a4':
                self.do_string(key, 4)
            elif self.FLAG[key]['format'] == '10I8':
                self.do_string(key, 8)
            elif self.FLAG[key]['format'] == '5E16.8':
                self.do_string(key, 16)
            elif self.FLAG[key]['format'] == '1a80':
                self.do_string(key, 80)
            else:
                exit(f"\n\tUNDEFINED format in {self.FLAG[key]['format']}\n")

    def do_string(self, key: str, split: int) -> list:
        """ fixing the data lists with FORTRAN format: 20a4, 10I8, 5E16.8, 1a80 """
        data_list = self.FLAG[key]['data']
        data_list = [data_list[i:i+split].strip() for i in range(0, len(data_list)-split+1, split)]
        self.FLAG[key]['data'] = data_list
        del data_list
    
    def mk_format(self) -> None:
        """correcting the data format """
        for key in self.FLAG.keys():
            if self.FLAG[key]['format'] == '10I8':
                self.do_integer(key)
            elif self.FLAG[key]['format'] == '5E16.8':
                self.do_float(key)
            else:
                pass

    def do_integer(self, key: int) -> None:
        # fixing the integer format of 10I8
        self.FLAG[key]['data'] = [int(x) for x in self.FLAG[key]['data']]
    
    def do_float(self, key: int) -> None:
        # fixing the float format of 5E16.8
        self.FLAG[key]['data'] = [float(x) for x in self.FLAG[key]['data']]


class GETTOP:
    """
    break down the data information read from TOP file
    This class only extract data from TOPFILE which was read by TOP class.
    Extracting data (anlaysing!), atom-type, ... to be used for making LAMMPS input file will be 
    done by DOTOP class
    """

    def __init__(self, top: dict) -> None:
        self.top = top
        del top
        print("  Setting the attributs ...")

    def get_top(self) -> None:
        # set general attributes which has either NATOM or NRES data-length
        self.set_attributes()
        # get genral data about types: masses, charges, ...
        self.get_types()

    def set_attributes(self) -> None:
        # call functions to set attributes
        self.get_pointers()
        self.get_atom_name()
        self.get_charges()
        self.get_masses()
        self.get_atom_type()
        self.get_residue_label()
        self.get_residue_pointer()
        self.get_LJ_coeff()
        self.get_bond_coeff()
        self.get_pair_index()
        self.print_info()

    def get_pointers(self) -> None:
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
        NPHB    :    Number of distinct 10-12 hydrogen bond pair types 2 IFPERT Set to 1 
         if topology contains residue perturbation information.
        NBPER   :    Number of perturbed bonds
        NGPER   :    Number of perturbed angles
        NDPER   :    Number of perturbed torsions
        MBPER   :    Number of bonds in which both atoms are being perturbed
        MGPER   :    Number of angles in which all 3 atoms are being perturbed
        MDPER   :    Number of torsions in which all 4 atoms are being perturbed 3 IFBOX Flag
         indicating whether a periodic box is present. Values can be 0 (no box),
        1 (orthorhombic box) or 2 (truncated octahedro
        NMXRS   :    Number of atoms in the largest residue IFCAP Set to 1 if a solvent CAP is being used
        NUMEXTRA:    Number of extra points in the topology file
        NCOPY   :    Number of PIMD slices or number of beads
        https://ambermd.org/prmtop.pdf
        """
        length = len(self.top['POINTERS']['data'])
        # setting the attributes:
        # set them to None
        nones = lambda n: [None for _ in range(n)]
        self.NATOM, self.NTYPES, self.NBONH, self.MBONA, self.NTHETH,\
            self.MTHETA, self.NPHIH, self.MPHIA, self.NHPARM, self.NPARM,\
            self.NNB, self.NRES, self.NBONA, self.NTHETA, self.NPHIA, self.NUMBND,\
            self.NUMANG, self.NPTRA, self.NATYP, self.NPHB, self.IFPERT, self.NBPER,\
            self.NGPER, self.NDPER, self.MBPER, self.MGPER, self.MDPER,\
            self.IFBOX, self.NMXRS, self.IFCAP, self.NUMEXTRA\
            = nones(31)
        # NCOPY may not be present!
        self.NCOPY = nones(1)
        # setting all the data
        self.NATOM, self.NTYPES, self.NBONH, self.MBONA, self.NTHETH,\
            self.MTHETA, self.NPHIH, self.MPHIA, self.NHPARM, self.NPARM,\
            self.NNB, self.NRES, self.NBONA, self.NTHETA, self.NPHIA, self.NUMBND,\
            self.NUMANG, self.NPTRA, self.NATYP, self.NPHB, self.IFPERT, self.NBPER,\
            self.NGPER, self.NDPER, self.MBPER, self.MGPER, self.MDPER,\
            self.IFBOX, self.NMXRS, self.IFCAP, self.NUMEXTRA\
            = self.top['POINTERS']['data'][0:31]
        if length > 31:
            self.NCOPY = self.top['POINTERS']['data'][31]

    def get_atom_name(self) -> None:
        """
        make a list of the atoms' name
        his section contains the atom name for every atom in the prmtop.
        %FORMAT(20a4) There are NATOM 4-character strings in this section.
        """
        length = len(self.top['ATOM_NAME']['data'])
        if length != self.NATOM:
            exit(f"NATOM != N of ATOM_NAME: {length}")
        self.ATOM_NAME = self.top['ATOM_NAME']['data']
        
    def get_charges(self) -> None:
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
        if length != self.NATOM:
            exit(f"NATOM != N of CHARGE: {length}")
        # convert to [e] unit
        charges = [q/kele for q in charges]
        # correct the precision
        charges = [np.round(q, decimals=10) for q in charges]
        self.top['CHARGE']['data'] = charges
        del charges
        
    def get_masses(self) -> None:
        """
        This section contains the atomic mass of every atom in g mol^{-1}.
        %FORMAT(5E16.8)
        There are NATOM floating-point numbers in this section.
        """
        masses = self.top['MASS']['data']
        length = len(masses)
        if length != self.NATOM:
            exit(f"NATOM != N of MASS: {length}")
        self.top['MASS']['data'] = masses
        del masses
    
    def get_atom_type(self) -> None:
        """
        This section contains the Lennard-Jones atom type index. The Lennard-
        Jones potential contains parameters for every pair of atoms in the system.
        %FORMAT(10I8)
        There are NATOM integers in this section.
        """
        atom_type = self.top['ATOM_TYPE_INDEX']['data']
        length = len(atom_type)
        if length != self.NATOM:
            exit(f"NATOM != N of ATOM_TYPE_INDEX: {length}")
        self.top['ATOM_TYPE_INDEX']['data'] = atom_type
        self.ATOM_TYPE_INDEX = atom_type
        del atom_type
    
    def get_residue_label(self) -> None:
        """
        This section contains the residue name for every residue in
        the prmtop. Residue names are limited to 4 letters, and might
        not be whitespace-delimited if any residues have 4-letter names.
        %FORMAT(20a4)
        There are NRES 4-character strings in this section
        """
        length = len(self.top['RESIDUE_LABEL']['data'])
        if length != self.NRES:
            exit(f"NRES != N of RESIDUE_LABEL: {length}")

    def get_residue_pointer(self) -> None:
        """
        This section lists the first atom in each residue.
        %FORMAT(10i8)
        There are NRES integers in this section
        """
        residue = self.top['RESIDUE_POINTER']['data']
        length = len(residue)
        if length != self.NRES:
            exit(f"NRES != N of RESIDUE_POINTER: {length}")
        self.top['RESIDUE_POINTER']['data'] = residue
        del residue

    def get_LJ_coeff(self) -> None:
        """
        This section contains the LJ A and B-coefficients (a_{i,j}, 
        b_{i,j} in Eq. LENNARD JONES) for all pairs of distinct LJ
        types (see sections ATOM TYPE INDEX and NONBONDED PARM INDEX).
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
        if len_lj != alength:
            exit("\n\tWRONG N of LJ A coeffs\n")
        if len_lj != blength:
            exit("\n\tWRONG N of LJ B coeffs\n")
        self.LJA = acoeffs
        self.LJB = bcoeffs
        self.top['LENNARD_JONES_ACOEF']['data'] = acoeffs
        self.top['LENNARD_JONES_BCOEF']['data'] = bcoeffs

    def get_types(self) -> None:
        """
        making a DataFrame from the atoms' name, charges, residues
        to extract the information about the number of types,
        with their symbols and properties
        """
        data_dict = dict()
        for key in self.top.keys():
            if len(self.top[key]['data']) == self.NATOM:
                data_dict[key] = self.top[key]['data']
        # make a dataframe to extract infos
        self.df = self.mk_df(data_dict)
        del data_dict

    def mk_df(self, data_dict) -> pd.DataFrame:
        df = pd.DataFrame.from_dict(data_dict)
        df = df.drop(['EXCLUDED_ATOMS_LIST'], axis=1)
        df.to_csv('df', sep='\t', index=False)
        df = df.groupby('ATOM_NAME', as_index=False).mean()
        df = df.astype({"ATOM_TYPE_INDEX": int})
        return df

    def get_bond_coeff(self) -> None:
        """
        Bond's coeffs for the harmonic interactions
        distance and K 
        E_{bond} = 0.5 * K * (r - r_{eq})^2
        """
        self.get_bond_force()
        self.get_bond_r_eq()

    def get_bond_force(self) -> None:
        """
        Bond energies are calculated according to the equation for the harmonic interaction
        This section lists all of the bond force constants (k in Eq. 2) in units
        kcal mol^{-1}  ̊A^{-2} for each unique bond type. Each bond in BONDS INC HYDROGEN
        and BONDS WITHOUT HYDROGEN (see below) contains an index into this array.
        %FORMAT(5E16.8)
        There are NUMBND floating point numbers in this section.
        """
        self.BOND_FORCE_CONSTANT= self.top['BOND_FORCE_CONSTANT']['data'] 
    
    def get_bond_r_eq(self) -> None:
        """
        This section lists all of the bond equilibrium distances (~req ) in
        units of Angestrom for each unique bond type. This list is indexed the same way as
        BOND FORCE CONSTANT.
        %FORMAT(5E16.8)
        There are NUMBND floating point numbers in this section
        """
        self.BOND_EQUIL_VALUE = self.top['BOND_EQUIL_VALUE']['data']

    def get_pair_index(self) -> None:
        """
        This section contains the pointers for each pair of LJ atom types into the
        LENNARD JONES ACOEF and LENNARD JONES BCOEF arrays (see below). The
        pointer for an atom pair in this array is calculated from the LJ atom type
        index of the two atoms (see ATOM TYPE INDEX above).
        %FORMAT(10I8)
        There are NTYPES * NTYPES integers in this section.
        """
        self.NONBONDED_PARM_INDEX = self.top['NONBONDED_PARM_INDEX']['data']

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
        with open(PDBFILE, 'r') as f:
            while True:
                line = f.readline()
                lineCounter += 1
                if line.strip().startswith('ATOM'):
                    line = line.strip().split(' ')
                    line = [item for item in line if item]
                    # ATOM      9 Si   Si      9      -6.272  -1.062  -5.117  1.00  0.00
                    _, i_id, i_name, i_residue, i_chain, i_x, i_y, i_z, _, _ = line
                    id.append(i_id)
                    name.append(i_name)
                    residue.append(i_residue)
                    chain.append(i_chain)
                    x.append(i_x)
                    y.append(i_y)
                    z.append(i_z)
                if not line:
                    break

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

    def mk_lmp_df(self,
                  id: int,
                  chain: int,
                  name: str,
                  q: float,
                  x: float,
                  y: float,
                  z: float,
                  nx: int,
                  ny: int,
                  nz: int,
                  sharp: str) -> dict:
        """
        making a DataFrame in LAMMPS 'full' atom style:

        id molecule-tag atom-type q x y z nx ny nz

        Since the atom type is not defined yet, the atom name will be used, and later it will be replaced with the atom type from TOPFILE. 
        Also, all q=0.0 will be set later the charges information will write into a different file.
        'chain' here is the same as molecule-tag in LAMMPS
        """
        data_dict = {'id': id,
                     'chain': chain,
                     'name': name,
                     'q': q,
                     'x': x,
                     'y': y,
                     'z': z,
                     'nx': nx,
                     'ny': ny,
                     'nz': nz,
                     'sharp': sharp,
                     'symbol': name}
        self.lmp_df = pd.DataFrame(data_dict)
        self.NATOM = len(id)
        self.NRES = np.max(chain)
        self.ALL_NAME = name
        self.ATOM_NAMES= list(set(name))
        self.NNAMES = len(self.ATOM_NAMES)
        del id, chain, name, q, x, y, z, nx, ny, nz, sharp, data_dict
        self.print_info()

    def print_info(self) -> typing.TextIO:
        print(f"\t seeing {self.NATOM}\t atoms")
        print(f"\t seeing {self.NRES}\t reseidues (molecules)")
        print(f"\t seeing {self.NNAMES}\t atom names")
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

    def __init__(self, top: dict, pdb: dict) -> None:
        self.top = top
        self.pdb = pdb
        del top

    def mk_bonds(self) -> None:
        # print(f"BOND_FORCE_CONSTANT: {self.top.top['BOND_FORCE_CONSTANT']['data']}")
        # print(f"BOND_EQUIL_VALUE: {self.top.top['BOND_EQUIL_VALUE']['data']}")
        # print(f"BONDS_WITHOUT_HYDROGEN: {self.top.top['BONDS_WITHOUT_HYDROGEN']['data']}")
        h_bonds = self.get_h_bond()
        self.bond_df = self.mk_hbond_df(h_bonds)
        self.set_attributes()

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
        h_bonds = [h_bonds[i: i+3] for i in range(0, len(h_bonds), 3)]
        h_bonds = [self.crct_index(lst) for lst in h_bonds]
        h_bonds = [self.append_name(lst) for lst in h_bonds]
        return h_bonds

    def crct_index(self, lst: list) -> list:
        for i in range(2):
            lst[i] = self.return_index(lst[i])
        return lst

    def return_index(self, x: int) -> int: return int((x/3)+1)

    def mk_hbond_df(self, h_bonds: list) -> pd.DataFrame:
        """return datafram from h_bond list"""
        columns = ['ai', 'aj', 'type', 'cmt', 'i_name', 'j_name']
        return pd.DataFrame(h_bonds, columns=columns)

    def set_attributes(self) -> None:
        # set attributes for bonnds
        self.NBONDS = len(self.bond_df)
        self.NBTYPES = self.bond_df['type'].max()
        self.bond_df.index += 1

    def append_name(self, lst: list) -> list:
        lst.append("#")
        lst.append(self.pdb.ALL_NAME[lst[0]-1])
        lst.append(self.pdb.ALL_NAME[lst[1]-1])
        return lst


class LMPPAIR:
    """
    Extract pair interaction between types (None bonded interactions)
    related cards to work with:
        - LENNARD JONES ACOEF
        - LENNARD JONES BCOEF
        - NONBONDED PARM INDEX
        - EXCLUDED ATOMS LIST
    """
    def __init__(self, top: dict, lmp: dict) -> None:
        self.top = top
        self.lmp = lmp
        del top, lmp

    def set_pairs(self) -> None:
        self.get_coeff()
        lst = self.get_index()
        self.Pair_df = self.mk_df(lst)
    
    def get_index(self) -> list:
        """
        Get indices for nonbonded interactions
        The index for two atoms i and j into the LENNARD JONES ACOEF and
        LENNARD JONES BCOEF arrays is calculated as:

        index = NONBONDED PARM INDEX [NTYPES * (ATOM TYPE INDEX(i) - 1) + ATOM TYPE INDEX(j)]

        Note, each atom pair can interact with either the standard 12-6 LJ po-
        tential or via a 12-10 hydrogen bond potential. If index in the equation is negative,
        then it is an index into HBOND ACOEF and HBOND BCOEF instead

        %FORMAT(10I8)
        There are NTYPES * NTYPES integers in this section.
        """
        lst = []
        # swap name and type in self.lmp.types
        type_name = {int(v): k for k, v in self.lmp.types.items()}
        pair_lst = [i for i in range(1, self.top.NTYPES+1)]
        for i, j in itertools.combinations_with_replacement(pair_lst, 2):
            try:
                inx = self.top.NTYPES * (i-1) + j
                ind = self.top.NONBONDED_PARM_INDEX[inx-1]
                # since Na+ is type 5, we drop the interaction here!
                # also append LAMMPS scripts for the pair interaction
                if i == 5 or j==5: 
                    pass
                else:
                    i_name = re.sub("\d", "", type_name[i])
                    j_name = re.sub("\d", "", type_name[j])
                    lst.append([i, j, self.sigma[ind-1], self.epsilon[ind-1],
                                'pair_coeff', 'lj/cut', '#', i_name, j_name])
            except:
                exit(f"\t WRONG interactions! between: {i} and {j}")
        return lst
    
    def get_coeff(self) -> None:
        self.sigma, self.epsilon = self.convert_unit()

    def convert_unit(self) -> list:
        """
        Converting A and B to epsilon and sigma
        AMBER  : E_{ij} = [a_{ij}/r^{12}] - [b_{ij}/r^{6}]
        LAMMPS : E_{ij} = 4 * epsilon [ (sigma/r)^{12} - (sigma/r)^{6} ]
        => sigma = (a / b)^(1 / 6)
           epsilon = b^2 / 4a
        """
        epsilon = []
        sigma = []
        for a, b in zip(self.top.LJA, self.top.LJB):
            sigma.append((a/b)**(1/6))
            epsilon.append(b**2 / (4*a))
        return sigma, epsilon

    def mk_df(self, lst) -> pd.DataFrame:
        columns=['ai', 'aj', 'sigma', 'epsilon', 'lmp_rule', 'style', '#', 'i_name', 'j_name']
        df = pd.DataFrame(lst, columns=columns)
        return df

    
class LMPDATA:
    """
    Update DataFrame for LAMMPS input:
    data file in full Atom style:
    id mol type q x y z nx ny nz
    """
    def __init__(self, pdb: dict, top: dict, bonds: dict) -> None:
        self.lmp_df = pdb.lmp_df
        self.top = top
        self.bonds = bonds
        del pdb, top, bonds

    def mk_lmp(self) -> None:
        self.types, self.charges, self.masses = self.get_q_name_mass()
        # there are Na+ in files, remove them and update all the related attributs
        self.data = self.update_df()
        self.ppo_key('Na+')
        # get box
        self.get_box()
        # write LAMMPS data file (DATAFILE)
        self.write_data()

    def ppo_key(self, popkey: str) -> None:
        # checking for the Na+
        try:
            # removing Na+
            self.lmp_df = self.lmp_df[self.lmp_df.symbol != popkey]
            self.NATOM = len(self.lmp_df.id)
            self.NTYPES = self.top.NTYPES - 1
            self.NRES = self.lmp_df.chain.max()
            # remove from masses, charges, and types
            self.types.pop(popkey)
            self.charges.pop(popkey)
            self.masses.pop(popkey)
        except:
            print(f"\tthere were no {popkey} to drop!!\n")

    def to_orgin(self) -> None:
        # put mins to origin (0,0,0)
        self.lmp_df.x = self.move_to_zero(self.lmp_df.x)
        self.lmp_df.y = self.move_to_zero(self.lmp_df.y)
        self.lmp_df.z = self.move_to_zero(self.lmp_df.z)

    # move the minumus to zero
    def move_to_zero(self, data: list) -> list: return data-np.min(data)

    def update_df(self) -> pd.DataFrame:
        """
        Update DataFram from PDB (self.lmp_df) by substituting the type and charge with data
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
        self.xlo = self.lmp_df['x'].min() - OFFSET
        self.xhi = self.lmp_df['x'].max() + OFFSET
        self.ylo = self.lmp_df['y'].min() - OFFSET
        self.yhi = self.lmp_df['y'].max() + OFFSET
        self.zlo = self.lmp_df['z'].min() - OFFSET
        self.zhi = self.lmp_df['z'].max() + OFFSET

    def write_data(self) -> typing.TextIO:
        with open(DATAFILE, 'w') as f:
            f.write(f"# dtat from: {TOPFILE} and {PDBFILE}\n")
            f.write(f"\n")
            f.write(f"{self.NATOM} atoms\n")
            f.write(f"{self.NTYPES} atom types\n")
            f.write(f"{self.bonds.NBONDS} bonds\n")
            f.write(f"{self.bonds.NBTYPES} bond types\n")
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
            columns = ['type', 'ai', 'aj', 'cmt', 'i_name', 'j_name']
            self.bonds.bond_df.to_csv(f, sep='\t', index=True, header=None,
                                      columns=columns, float_format='%g')
            self.print_info()

    def print_info(self) -> typing.TextIO:
        print(f"Writing '{DATAFILE}' (LAMMPS data file) ... \n")
        print(f"\t writting {self.NATOM}\t atoms")
        print(f"\t writting {self.NRES}\t reseidues (molecules)")
        print(f"\t writting {self.NTYPES}\t atom types")
        print(f"\t writting {self.bonds.NBONDS}\t bonds")
        print(f"\t writting {self.bonds.NBTYPES}\t bond types")
        print(f"\n")


class LMPPARAM:
    """
    write out parameters about the system such as
    mass, coeffs, bonds, ...
    """
    def __init__(self, lmp: dict, bond: dict, pair: dict) -> None:
        self.lmp = lmp
        self.bond = bond
        self.pair = pair
        del lmp, pair

    def mk_types(self) -> None:
        # call functions one by one to set attributes
        self.get_types()
        self.get_masses()
        self.get_charges()
        self.write_parameters()
        self.print_info()

    def get_types(self) -> None:
        """
        make clean set of the types, becuase NTYPES != NNAMES
        """
        # convert the the format of the types: int -> str
        self.lmp.types = {k: str(v) for k, v in self.lmp.types.items()}
        # swap keys, and values them tuples the new key
        d = {tuple(v): k for k, v in self.lmp.types.items()}
        # swap the key and value again
        self.lmp.types = {v: k[0] for k, v in d.items()}
        del d

    def get_masses(self) -> None:
        """ Update masses attribute based on the types """
        self.lmp.masses = {k: v for k, v in self.lmp.masses.items()
                           if k in self.lmp.types.keys()}

    def get_charges(self) -> None:
        """ Update charges (NOT needed, bu incase!) attribute based
        on the types
        """
        self.lmp.charges = {k: v for k, v in self.lmp.charges.items()
                            if k in self.lmp.types.keys()}

    def write_parameters(self) -> typing.TextIO:
        with open(PARAMFILE, 'w') as f:
            f.write(f"# Parameters for '{DATAFILE}' \
                from '{TOPFILE}' and '{PDBFILE}'\n")
            f.write(f"\n")
            self.write_mass(f)
            self.write_q(f)
            self.write_group(f)
            self.write_pair(f)
            self.write_bond(f)
            self.write_constrains(f)

    def write_mass(self, f: typing.IO) -> typing.TextIO:
        f.write(f"# mass of each type\n")
        for key, value in self.lmp.masses.items():
            f.write(f"mass {self.lmp.types[key]}\t{value:.3f}\t#\
                     {self.drop_digit(key)}\n")
        f.write("\n")

    def write_q(self, f: typing.IO) -> typing.TextIO:
        f.write(f"# charge of each type\n")
        f.write(f"# charges are already set in data file ({DATAFILE}),\
                 here added as comments\n")
        f.write(f"CC\n")
        for key, value in self.lmp.charges.items():
            f.write(f"set type {self.lmp.types[key]} \
                charge {value:.3f}\t# {self.drop_digit(key)}\n")
        f.write(f"CC\n")
        f.write(f"\n")

    def write_group(self, f: typing.IO) -> typing.TextIO:
        # define Hydrogen group
        f.write(f"# define name of the group \n")
        f.write(f"group Hydrogen type {self.lmp.types['H']}\n")
        f.write(f"group OH type {self.lmp.types['OH']}\n")
        f.write(f"group OM type {self.lmp.types['OM3']}\n")
        f.write(f"group Silicon type {self.lmp.types['Si']}\n")
        f.write(f"group Oxygen type {self.lmp.types['OH']}\
            {self.lmp.types['OM3']}\n")
        f.write(f"group freeze type {self.lmp.types['OH']}\
            {self.lmp.types['OM3']} {self.lmp.types['Si']}\n")
        f.write(f"\n")

    def write_pair(self, f: typing.IO) -> typing.TextIO:
        # write ij pair interactions
        columns = ['lmp_rule', 'ai', 'aj', 'style',
                   'epsilon', 'sigma', '#', 'i_name', 'j_name']
        f.write(f"# define the interactions between particles\n")
        f.write(f"# {'   '.join(columns)}\n")
        self.pair.Pair_df.to_csv(f, sep='\t', columns=columns,
                                 header=None, index=False)
        f.write(f"\n")

    def write_bond(self, f: typing.IO) -> typing.TextIO:
        # writting bond coeffs
        f.write(f"# bond coeff for the all the atom types \n")
        for i in range(self.bond.NBTYPES):
            f.write(f"bond_coeff {i+1} harmonic")
            f.write(f"{self.lmp.top.BOND_FORCE_CONSTANT[i]} {self.lmp.top.BOND_EQUIL_VALUE[i]}\n")
        f.write(f"\n")

    def write_constrains(self, f) -> typing.TextIO:
        f.write(f"# set zero forces on the slab group beside hydrogen\n")
        f.write(f"fix freeze freeze setforce 0.0 0.0 0.0\n")
        f.write(f"\n")

    def print_info(self) -> typing.TextIO:
        print(f"Writing parameters in '{PARAMFILE}' ...\n")

    def drop_digit(self, obj: str) -> str: return re.sub("\d", "", obj)


if __name__ == "__main__":
    TOPFILE = "test3.top"
    PDBFILE = "test.pdb"
    DATAFILE = "slab.data"
    PARAMFILE = 'parameters.lmp'
    OFFSET = 1
    data = READTOP()
    data.get_data()
    top = GETTOP(data.FLAG)
    top.get_top()
    pdb = PDB()
    pdb.read_pdb()
    lmpbond = LMPBOND(top, pdb)
    lmpbond.mk_bonds()
    lmpdata = LMPDATA(pdb, top, lmpbond)
    lmpdata.mk_lmp()
    lmppair = LMPPAIR(top, lmpdata)
    lmppair.set_pairs()
    lmpparam = LMPPARAM(lmpdata, lmpbond, lmppair)
    lmpparam.mk_types()
