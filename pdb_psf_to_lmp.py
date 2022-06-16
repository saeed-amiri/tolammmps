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
                    data_list.append(self.process_line(line))
                if not line:
                    break
        return data_list

    def process_line(self, line: str) -> list:
        """Process the line on based on the PDB file"""
        # first check the length of the line MUST be equal to 79
        if len(line) != 79:
            exit(f"ERROR! wrong line length: {len(line)} != 79")
        atom_id: int = int(line[6:11].strip())
        atom_name: str = line[13:16].strip()
        residue_name: str = line[17:20].strip()
        residue_number: str = line[22:27].strip()
        atom_x: float = float(line[30:39].strip())
        atom_y: float = float(line[39:47].strip())
        atom_z: float = float(line[47:55].strip())
        atom_symbool: str = line[76:78].strip()
        return [atom_id,
                atom_name,
                residue_name,
                residue_number,
                atom_x,
                atom_y,
                atom_z,
                atom_symbool]

    def mk_df(self, data: list[list]) -> pd.DataFrame:
        """Making DataFrame from PDB file"""
        columns: list[str] = ['atom_id',
                              'atom_name',
                              'residue_name',
                              'residue_number',
                              'atom_x',
                              'atom_y',
                              'atom_z',
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

    This file contains information about bonds, angles, dihedrals, ...

    """


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
