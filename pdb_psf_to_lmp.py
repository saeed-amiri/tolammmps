import sys
import typing
import pandas as pd
import numpy as np


class Doc:
    """ Convert PDB  and PSF file to LAMMPS data file
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


class PDB:
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
        self.read_pdb()

    def read_pdb(self) -> None:
        """Reading line by line of the pdb file"""
        with open(PDBFILE, 'r') as f:
            while True:
                line = f.readline()
                if line.strip().startswith("ATOM"):
                    self.process_line(line)
                if not line:
                    break

    def process_line(self, line) -> None:
        """ Process the line on based on the PDB file"""
        # first check the length of the line MUST be equal to 79
        if len(line) != 79:
            exit(f"ERROR! wrong line length: {len(line)} != 79")
        atom_id: int = int(line[6:11])
        atom_name: str = line[13:16]
        residue_name: str = line[17:20]
        residue_number: str = line[22:27]
        atom_x: float = float(line[30:39])
        atom_y: float = float(line[39:47])
        atom_z: float = float(line[47:55])
        atom_symbool: str = line[76:78]



# Constant values
if len(sys.argv) < 2:
    doc = Doc()
    exit(doc.__doc__)
    
PDBFILE = sys.argv[1].__add__('.pdb')
PSFFILE = sys.argv[1].__add__('.psf')
LMPFILE = sys.argv[1].__add__('.data')

pdb = PDB()
pdb.get_data()
