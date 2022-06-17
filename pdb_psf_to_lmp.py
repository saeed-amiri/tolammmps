from pprint import pprint
import re
import sys
import typing
import pandas as pd
import numpy as np


class Doc:
    """Convert PDB  and PSF file to LAMMPS data file
    Input:
        Name of the file without extension (sys.argv[1])
        Both PDB and PSF should have the same name
    Output:
        LAMMPS data file
        [ checking PEP8
        ~/.local/bin/pycodestyle pdb_psf_to_lmp.py
        and typing"
        ~/.local/bin/mypy pdb_psf_to_lmp.py
    ]

    USAGE example: pdb_psf_to_lmp.py solavte

    """


class Pdb:
    """
    Reading the PDB file and returning the coordinates.
    The PDB file consider is in standard PDB format.
    It reads the PDB file and return the information in DataFrames.
    convert LAMMPS data file to a standard PDB file format based on:
    [https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html]
    A PDB file has a format as:
    ATOM 	atomic coordinate record containing the X,Y,Z
        orthogonal Å coordinates for atoms in standard residues]
        (amino acids and nucleic acids).
    Protein Data Bank Format:
      Coordinate Section
        Record Type	Columns	Data 	Justification	Data Type
        ATOM 	1-4	“ATOM”	                    	character
        7-11#	Atom serial number      	right	integer
        13-16	Atom name	                left*	character
        17	    Alternate location indicator		character
        18-20§	Residue name	            right	character
        22	    Chain identifier		            character
        23-26	Residue sequence number	    right	integer
        27	    Code for insertions of residues		character
        31-38	X orthogonal Å coordinate	right	real (8.3)
        39-46	Y orthogonal Å coordinate	right	real (8.3)
        47-54	Z orthogonal Å coordinate	right	real (8.3)
        55-60	Occupancy	                right	real (6.2)
        61-66	Temperature factor	        right	real (6.2)
        73-76	Segment identifier¶	        left	character
        77-78	Element symbol              right	character
        79-80	Charge		                        character

    #Chimera allows (nonstandard) use of columns 6-11 for the integer
        atom serial number in ATOM records, and in TER records, only the
        “TER” is required.

    *Atom names start with element symbols right-justified in columns
        13-14 as permitted by the length of the name. For example, the
        symbol FE for iron appears in columns 13-14, whereas the symbol
        C for carbon appears in column 14 (see Misaligned Atom Names).
        If an atom name has four characters, however, it must start in
        column 13 even if the element symbol is a single character
        (for example, see Hydrogen Atoms).

    §Chimera allows (nonstandard) use of four-character residue names
        occupying an additional column to the right.

    ¶Segment identifier is obsolete, but still used by some programs.
        Chimera assigns it as the atom attribute pdbSegment to allow
        command-line specification.

    The format of ecah section is (fortran style):
    Format (A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)
    """

    def __init__(self) -> None:
        print(f"Reading '{PDBFILE}' ...")

    def get_data(self) -> None:
        """Read and get the atomic strcuture in lammps`"""
        # Read PDB file and return a list of list
        data = self.read_pdb()
        # Convert the list to DataFrame
        df = self.mk_df(data)
        # Check the residue number order
        self.atoms_df = self.check_residue_number(df)

    def read_pdb(self) -> list:
        """Reading line by line of the pdb file"""
        data_list: list[typing.Any] = []
        with open(PDBFILE, 'r') as f:
            while True:
                line = f.readline()
                if line.strip().startswith("ATOM"):
                    data_list.append(self.__process_line(line))
                if not line:
                    break
        return data_list

    def __process_line(self, line: str) -> list:
        """Process the line on based on the PDB file"""
        # first check the length of the line MUST be equal to 79
        if len(line) != 79:
            exit(f"ERROR! wrong line length: {len(line)} != 79")
        atom_id: int = int(line[6:11].strip())
        atom_name: str = line[13:16].strip()
        residue_name: str = line[17:20].strip()
        residue_number: str = line[22:27].strip()
        x: float = float(line[30:39].strip())
        y: float = float(line[39:47].strip())
        z: float = float(line[47:55].strip())
        atom_symbool: str = line[76:78].strip()
        return [atom_id,
                atom_name,
                residue_name,
                residue_number,
                x,
                y,
                z,
                atom_symbool]

    def mk_df(self, data: list[list]) -> pd.DataFrame:
        """Making DataFrame from PDB file"""
        columns: list[str] = ['atom_id',
                              'atom_name',
                              'residue_name',
                              'residue_number',
                              'x',
                              'y',
                              'z',
                              'atom_symbool']
        df = pd.DataFrame(data, columns=columns)
        return df

    def check_residue_number(self, df: pd.DataFrame) -> pd.DataFrame:
        """The residue number in the PDB file created by VMD usually
        does not have the correct order, and the numbering is almost
        random.
        Here we check that and if needed will be replaced with correct
        one.
        """
        # Make sure the residue type is integer
        df = df.astype({'residue_number': int})
        # Get the uniqe numbers of residues
        df_reduced: pd.DataFrame = df.groupby(df.residue_number).mean()
        # add a new index from 1
        df_reduced = df_reduced.reset_index()
        df_reduced.index += 1
        # make a dict from the new index and the residue
        residue_dict: dict[int, int] = {
            k: v for k, v in zip(df_reduced.residue_number, df_reduced.index)
            }
        # Make a list of new residue number for the main DataFrame\
        # I call it 'mol'
        mol: list[int] = []
        for _, row in df.iterrows():
            mol.append(residue_dict[row['residue_number']])
        # Add the new mol numbers to the DataFrame
        df['mol'] = mol
        del df_reduced
        del residue_dict
        return df


class Psf:
    """Read PSF file
    www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node23.html
    A PSF file, also called a protein structure file, contains all of
    the molecule-specific information needed to apply a particular force
    field to a molecular system.
    The PSF file contains six main sections of interest:
        atoms, bonds, angles, dihedrals, impropers and cross-terms.
    Each section starts with the number of data in that section and
    name of the section which begins with: !.
    This file contains information about bonds, angles, dihedrals, ...
    """

    def __init__(self) -> None:
        pass

    def get_data(self) -> None:
        """Since the number of sections is fixed, here, the class
        explicitly read the data instead of figuring out where and how
        much info there is."""
        attr_list = self.get_sections()
        self.read_data(attr_list)

    def get_sections(self) -> list[str]:
        """Since the name of the sections are fixed, set the attributs
        for number of data in each section
        """
        attr_list: list[str] = []
        with open(PSFFILE, 'r') as f:
            while True:
                line = f.readline()
                attr = self.set_section_attrs(line.strip())
                if attr is not None:
                    attr_list.append(attr)
                if not line:
                    break
        return attr_list

    def set_section_attrs(self, line: str) -> str:
        if "!" in line:
            l_line: list[str] = line.split(" ")
            l_line = [item for item in l_line if item]
            # Get the number of the each section
            if l_line[1].startswith("!"):
                # Keep the letters only get the name
                attr: str = re.sub('[^a-zA-Z]+', '', l_line[1])
                # Get the value
                value: int = int(l_line[0])
                # Set the attribute
                setattr(self, attr, value)
            else:
                attr = None
            return attr

    def read_data(self, attr_list: list[str]) -> None:
        """Reading the data"""
        # Making a dictionary of flages for the name of each section
        flages: dict[str: bool] = {k: False for k in attr_list}
        with open(PSFFILE, 'r') as f:
            while True:
                line = f.readline()
                flages = self.__process_line(line.strip(), flages)
                if not line:
                    break

    def __process_line(self,
                       line: str,
                       flages: dict[str: bool]) -> None:
        """Read the lines of the file and send each section to the
        proper function"""
        if "!" in line:
            l_line = line.split(' ')
            l_line = [item for item in l_line if item]
            attr: str = re.sub('[^a-zA-Z]+', '', l_line[1].strip())
            # Set all the flages to False and check for empty strings
            flages = {k: False for k in flages if k}
            # Set the correspond flag to True
            flages[attr] = True
        else:
            for k in flages:
                if k and flages[k]:
                    # Make a string to call the correspond function
                    # Drop the initit N
                    func_name = k[1:]
                    func_name = func_name.lower()
                    func_name = ''.join(['get_', func_name])
                    print(func_name, end=',')
        return flages

    def get_atom(self, line: str) -> None:
        pass

class WriteLmp:
    """Write the data in a full atoms style for LAMMPS
    Input:
        atoms_df (DataFrame from PDBFILE: Pdb class)
        bonds_df, angles_df, dihedrals, impropers_df (DataFrame from
        PSFFILE Psf class)
    Output:
        A LAMMPS data file
    """

    def __init__(self, atoms: pd.DataFrame) -> None:
        self.atoms = atoms
        del atoms

    def mk_lmp(self) -> None:
        """calling function to write data into a file"""
        # find box sizes
        self.set_box()
        # get number of atoms, types, bonds
        # self.set_numbers()
        # write file
        # self.write_data()
        # print(self.atoms['charge'].sum())

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
        with open(LMPFILE, 'w') as f:
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


if len(sys.argv) < 2:
    doc = Doc()
    exit(doc.__doc__)

# Constant values
PDBFILE = sys.argv[1].__add__('.pdb')
PSFFILE = sys.argv[1].__add__('.psf')
LMPFILE = sys.argv[1].__add__('.data')

if __name__ == '__main__':
    pdb = Pdb()
    pdb.get_data()
    psf = Psf()
    psf.get_data()
    print(psf.NATOM)
    lmp = WriteLmp(pdb.atoms_df)
    lmp.mk_lmp()
