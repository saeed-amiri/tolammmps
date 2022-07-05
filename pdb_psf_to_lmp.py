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
    It gets data implicitly, it reads the file and then based the data
    format it set the informations.
    For example the constants like : "NBOND", "NTHETA", "NATOM", "NBOND",
    and "Masses" it not get directly.
    It is just an experience :))

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
        atom_symbol: str = line[76:78].strip()
        return [atom_id,
                atom_name,
                residue_name,
                residue_number,
                x,
                y,
                z,
                atom_symbol]

    def mk_df(self, data: list[list]) -> pd.DataFrame:
        """Making DataFrame from PDB file"""
        columns: list[str] = ['atom_id',
                              'atom_name',
                              'residue_name',
                              'residue_number',
                              'x',
                              'y',
                              'z',
                              'atom_symbol']
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
        print(f"Reading '{PSFFILE}' ...")
        self.get_data()

    def get_data(self) -> None:
        """Since the number of sections is fixed, here, the class
        explicitly read the data instead of figuring out where and how
        much info there is."""
        attr_list = self.get_sections()
        self.read_data(attr_list)
        # For sake of simplicty I do this here!:
        self.attr_list = attr_list

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
        # Initiat a dictionary of flages for the name of each section
        flages: dict[str, bool] = {k: False for k in attr_list if k}
        # Initiat a dictionary of list to get data of each section
        self.data: dict[str, list[typing.Any]] =\
            {k: [] for k in attr_list if k}
        # Read PSFFILE line by line
        with open(PSFFILE, 'r') as f:
            while True:
                line = f.readline()
                if line.strip().startswith("REMARKS"):
                    pass
                else:
                    flages = self.__process_line(line.strip(), flages)
                if not line:
                    break

    def __process_line(self,
                       line: str,
                       flages: dict[str, bool]) -> dict[str, bool]:
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
                # check if key is not empty str, flage is True, and
                # line is not epmty
                if k and flages[k] and line.split():
                    self.data[k].append(line.split())
        return flages


class PsfToDf(Psf):
    """make DataFrame from data read by Psf
    Input:
        Dictionary of data from Psf
    Output:
        DataFrame from section for using in LAMMPS data file
    """

    def __init__(self) -> None:
        super().__init__()
        self.mk_dataframe()
        del self.data

    def mk_dataframe(self) -> None:
        for k in self.attr_list:
            # call the correspond function by name of the key
            func = getattr(self, self.mk_func_name(k))
            func(self.data[k])

    def mk_func_name(self, key: str) -> str:
        # Make a string to call the correspond function
        # Drop the initiate N
        func_name = key[1:]
        func_name = func_name.lower()
        func_name = ''.join(['get_', func_name])
        return func_name

    def get_title(self, data: list[str]) -> None:
        df = pd.DataFrame(data)
        if not df.empty:
            self.title = df

    def get_atom(self, data: list[str]) -> None:
        """The fields in the atom section are atom ID, segment name,
        residue ID, residue name, atom name, atom type, charge, mass,
        and an unused 0.

        Format of each line:
        int, str, int, str, str, str, float, float, int
        """
        dtype: dict[str, typing.Type[typing.Any]] =\
            {'atom_id': int,
             'segment_name': str,
             'residue_number': int,
             'residue_name': str,
             'atom_name': str,
             'atom_type': str,
             'charge': float,
             'mass': float,
             'unused': int}
        columns: list[str] = list(dtype.keys())
        df = pd.DataFrame(data, columns=columns)
        df = df.astype(dtype=dtype)
        df = self.check_residue_number(df)
        self.atoms_info = df
        del df

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

    def get_bond(self, data: list[str]) -> None:
        """The covalent bond section lists four pairs of atoms per line"""
        # first flatten the list
        bond_list: list[str] = [i for j in data for i in j]
        # An empty list for appending the lists in it
        bonds: list[list[str]] = []
        # There are NBOND bonds in the cards, which are pairs, so we
        # need to multiply it by 2 to get every individual atom id
        for i in range(0, 2*self.NBOND-1, 2):
            bonds.append([bond_list[i], bond_list[i+1]])
        columns = ['ai', 'aj']
        self.bonds = pd.DataFrame(bonds, columns=columns)
        del bonds

    def get_theta(self, data: list[str]) -> None:
        """The angle section lists three triples of atoms per line"""
        # first flatten the list
        angle_list: list[str] = [i for j in data for i in j]
        # An empty list for appending the lists in it
        angles: list[list[str]] = []
        # There are NBOND bonds in the cards, which are pairs, so we
        # need to multiply it by 3 to get every individual atom id
        for i in range(0, 3*self.NTHETA-1, 3):
            angles.append([angle_list[i], angle_list[i+1], angle_list[i+2]])
        columns = ['ai', 'aj', 'ak']
        self.angles = pd.DataFrame(angles, columns=columns)
        del angles

    def get_phi(self, data: list[str]) -> None:
        pass

    def get_imphi(self, data: list[str]) -> None:
        pass

    def get_don(self, data: list[str]) -> None:
        pass

    def get_acc(self, data: list[str]) -> None:
        pass

    def get_nb(self, data: list[str]) -> None:
        pass


class PsfToLmp(PsfToDf):
    """update all the data fram in the LAMMPS version
    Parent class:
        PsfToDf: contains all the DataFrame from PSF file
    Input:
        Atoms coordinate DataFrame from PDB
    Output:
        DataFrame to write data into LAMMPS file
    """

    def __init__(self, atoms: pd.DataFrame) -> None:
        super().__init__()
        self.atoms = atoms
        print(f"Converting to LAMMPS ...")
        del atoms

    def set_attrs(self) -> None:
        self.typ: pd.DataFrame = self.set_type()  # Include mass & q
        # Set the number of atom types
        self.Natom_type: int = len(self.typ)
        self.lmp_atoms: pd.DataFrame = self.mk_lmp_atoms()
        self.lmp_bonds: pd.DataFrame = self.mk_lmp_bonds()
        # Set the number of types
        self.Nbonds_type: int = self.lmp_bonds['typ'].max()
        self.lmp_angles: pd.DataFrame = self.mk_lmp_angles()
        self.Nangles_type: int = self.lmp_angles['typ'].max()

    def set_type(self) -> pd.DataFrame:
        """set the mass for each type"""
        columns: list[str] = ['charge', 'mass']
        df_sub: pd.DataFrame = self.atoms_info[columns].copy()
        # Copy the symbol of each atom to get a uniq column
        df_sub['atom_symbol'] = self.atoms['atom_symbol']
        df_sub = df_sub.groupby(by=['atom_symbol']).min()
        df_sub = df_sub.reset_index()
        df_sub.index += 1
        df_sub['typ'] = df_sub.index
        df_sub = df_sub.set_index('atom_symbol')
        return df_sub

    def mk_lmp_atoms(self) -> pd.DataFrame:
        """make atoms DataFrame in LAMMPS format"""
        # Columns to copy from PDB data for lammps
        columns: list[str] = ['atom_id', 'mol',
                              'x', 'y', 'z', 'atom_symbol']
        # Initiat the DataFrame
        df: pd.DataFrame = self.atoms[columns].copy()
        # Get the charges for each atom
        l_charge: list[float] = [
            self.typ.loc[item]['charge'] for item in self.atoms['atom_symbol']
        ]
        # Get the type for each atom
        l_typ: list[int] = [
            self.typ.loc[item]['typ'] for item in self.atoms['atom_symbol']
        ]
        l_typ = [int(typ) for typ in l_typ]
        df['charge'] = l_charge
        df['typ'] = l_typ
        # Add box flages
        df['nx'] = [0 for _ in l_typ]
        df['ny'] = [0 for _ in l_typ]
        df['nz'] = [0 for _ in l_typ]
        # Add # for comment information in the data file
        df['cmt'] = ['#' for _ in l_typ]
        # correct the index
        df.index += 1
        return df

    def mk_lmp_bonds(self) -> pd.DataFrame:
        """make DataFrame for LAMMPS Bonds section"""
        df: pd.DataFrame = self.bonds
        df = df.astype({'ai': int, 'aj': int})
        # Create the bonds for the atoms ai and aj
        # The index of lmp_atoms starts from 1 !!
        l_ai: list[str] = [
            self.lmp_atoms['atom_symbol'][item] for item in df['ai']
        ]
        l_aj: list[str] = [
            self.lmp_atoms['atom_symbol'][item] for item in df['aj']
        ]
        bond: list[str] = [f"{i}-{j}" for i, j in zip(l_ai, l_aj)]
        df['cmt'] = ['#' for _ in bond]
        df['bond'] = bond
        # Add id column
        df.index += 1
        df['id'] = df.index
        # Add type column
        df = self.get_bond_typ(df)
        return df

    def get_bond_typ(self, df: pd.DataFrame) -> pd.DataFrame:
        """Find how many differnt bonds are in the system"""
        columns: list[str] = ['bond']
        df_sub: pd.DataFrame = df[columns].copy()
        df_sub = df_sub.groupby(by=['bond']).min()
        df_sub = df_sub.reset_index()
        df_sub.index += 1
        df_sub['typ'] = df_sub.index
        df_sub = df_sub.set_index('bond')
        df['typ'] = [df_sub.loc[item]['typ'] for item in df['bond']]
        return df

    def mk_lmp_angles(self) -> pd.DataFrame:
        """make DataFrame for LAMMPS Angels section"""
        df: pd.DataFrame = self.angles
        df = df.astype({'ai': int, 'aj': int, 'ak': int})
        # Create the angles for the atoms ai and aj
        l_ai: list[str] = [
            self.lmp_atoms['atom_symbol'][item] for item in df['ai']
        ]
        l_aj: list[str] = [
            self.lmp_atoms['atom_symbol'][item] for item in df['aj']
        ]
        l_ak: list[str] = [
            self.lmp_atoms['atom_symbol'][item] for item in df['ak']
        ]
        angle: list[str] = [
            f"{i}-{j}-{k}" for i, j, k in zip(l_ai, l_aj, l_ak)
        ]
        df['cmt'] = ['#' for _ in angle]
        df['angle'] = angle
        columns: list[str] = ['angle']
        df_sub: pd.DataFrame = df[columns].copy()
        df_sub = df_sub.groupby(by=columns).min()
        df_sub = df_sub.reset_index()
        df_sub.index += 1
        df_sub['typ'] = df_sub.index
        df_sub = df_sub.set_index('angle')
        df['typ'] = [df_sub.loc[item]['typ'] for item in df['angle']]
        df.index += 1
        return df


class WriteLmp:
    """Write the data in a full atoms style for LAMMPS
    Input:
        atoms_df (DataFrame from PDBFILE: Pdb class)
        bonds_df, angles_df, dihedrals, impropers_df (DataFrame from
        PSFFILE Psf class)
    Output:
        A LAMMPS data file
    """

    def __init__(self, lmp: PsfToLmp) -> None:
        self.atoms = lmp.lmp_atoms
        self.bonds = lmp.lmp_bonds
        self.angles = lmp.lmp_angles
        self.Natoms = lmp.NATOM
        self.Natoms_type = lmp.Natom_type
        self.Nbonds = lmp.NBOND
        self.Nbonds_type = lmp.Nbonds_type
        self.Nangles = lmp.NTHETA
        self.Nangles_type = lmp.Nangles_type
        self.typ = lmp.typ
        print(f"Writting '{LMPFILE}' ...")

    def mk_lmp(self) -> None:
        """calling function to write data into a file"""
        # find box sizes
        self.set_box()
        # get number of atoms, types, bonds
        # self.set_numbers()
        # write file
        self.write_data()
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
            f.write(f"Data file from VDM for silica slab\n")
            f.write(f"\n")
            self.write_numbers(f)
            self.write_box(f)
            self.write_masses(f)
            self.write_atoms(f)
            self.write_bonds(f)
            self.write_angles(f)

    def write_numbers(self, f: typing.TextIO) -> None:
        f.write(f"{self.Natoms} atoms\n")
        f.write(f"{self.Natoms_type} atom types\n")
        f.write(f"{self.Nbonds} bonds\n")
        f.write(f"{self.Nbonds_type} bond types\n")
        f.write(f"{self.Nangles} angles\n")
        f.write(f"{self.Nangles_type} angle types\n")
        f.write(f"\n")

    def write_box(self, f: typing.TextIO) -> None:
        f.write(f"{self.xlim[0]:.3f} {self.xlim[1]:.3f} xlo xhi\n")
        f.write(f"{self.ylim[0]:.3f} {self.ylim[1]:.3f} ylo yhi\n")
        f.write(f"{self.zlim[0]:.3f} {self.zlim[1]:.3f} zlo zhi\n")
        f.write(f"\n")

    def write_masses(self, f: typing.TextIO) -> None:
        f.write(f"Masses\n")
        f.write(f"\n")
        for item in self.typ.iterrows():
            f.write(f"{int(item[1]['typ']): 3} {item[1]['mass']:.5f}\
                # {item[0]}\n")
        f.write(f"\n")

    def write_atoms(self, f: typing.TextIO) -> None:
        """Write atoms section inot file"""
        f.write(f"Atoms # full\n")
        f.write(f"\n")
        columns = ['mol', 'typ', 'charge', 'x', 'y', 'z', 'nx', 'ny', 'nz',
                   'cmt', 'atom_symbol']
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
        f.write(f"\n")

    def write_angles(self, f: typing.TextIO) -> None:
        f.write(f"Angles\n")
        f.write(f"\n")
        columns = ['typ', 'ai', 'aj', 'ak']
        self.angles.to_csv(f, sep=' ', index=True, columns=columns,
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
    psf = PsfToLmp(pdb.atoms_df)
    psf.set_attrs()
    lmp = WriteLmp(psf)
    lmp.mk_lmp()
