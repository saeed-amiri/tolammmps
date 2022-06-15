import typing
import pandas as pd
import numpy as np
import sys


class Doc:
    """ Convert PDB file to lammps data file
    Input:
        PDB file (sys.argv[1])
    Output:
        LAMMPS data file
        [ checking PEP8
        ~/.local/bin/pycodestyle pdb_to_lmp.py
        and typing"
        ~/.local/bin/mypy pdb_to_lmp.py
    ]
    """


class PDB:
    """
    reading the PDB file and returning the coordinates
    The PDB file is written in a normal format it is NOT in a standard
    PDB structure
    """
    def __init__(self) -> None:
        print(f"Reading '{PDBFILE}' ...")

    def get_data(self) -> None:
        """Read and get the atomic strcuture in lammps`"""
        pass
    def read_pdb(self) -> None:
        id, name, residue, chain, x, y, z = [], [], [], [], [], [], []
        lineCounter = 0
        with open(PDBFILE, 'r') as f:
            while True:
                line = f.readline()
                lineCounter += 1
                if line.strip().startswith('ATOM'):
                    line = line.strip().split(' ')
                    line = [item for item in line if item]
                    # ATOM   9 Si Si  9  -6.272 -1.062  -5.117  1.00  0.00
                    _, i_id, i_name, i_residue, i_chain, i_x, i_y, i_z, _, _ =\
                        line
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
                  sharp: str) -> None:
        """
        making a DataFrame in LAMMPS 'full' atom style:
        id molecule-tag atom-type q x y z nx ny nz
        Since the atom type is not defined yet, the atom name will be
        used, and later it will be replaced with the atom type from
        TOPFILE.
        Also, all q=0.0 will be set later the charges information will
        write into a different file.
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
        self.ATOM_NAMES = list(set(name))
        self.NNAMES = len(self.ATOM_NAMES)
        del id, chain, name, q, x, y, z, nx, ny, nz, sharp, data_dict
        self.print_info()

    def print_info(self) -> None:
        print(f"\t seeing {self.NATOM}\t atoms")
        print(f"\t seeing {self.NRES}\t reseidues (molecules)")
        print(f"\t seeing {self.NNAMES}\t atom names")
        print(f"\n")


PDBFILE = sys.argv[1]
pdb = PDB()