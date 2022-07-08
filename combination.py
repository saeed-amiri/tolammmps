import sys
import os
import re
import typing
import numpy as np
import pandas as pd


class Doc:
    """This code combines data files and returns a superstructure in
    LAMMPS full atom style.   The input file is <system>.struct file
    which contains the names and path of the  files and also a matrix
    that shows the order of the superstructure:

   ex. of an input file:

    First, it should have a symbol for each file and its path.
    The symbol is better to be an upper letter (for now, it does not
    recognize lower/upper case of the letter):

        ! D=decane.data
        ! W=water.data
        ! S=sio2.data
    Then it should have a matrix that shows how you want to build your
    superstructure:
        DWD
        _S_
        DWD
    which stack the blocks as: "decane water decane" most top level,
    SiO2 will stack in the second layer and "decane water decane" in
    the lowest layer. The atoms id will be set in the Z way of the
    matrix.

    Few reserved char:
    # (sharp): comment
    ! : to show the symbol and name of the files
    | (pipe char): which means the upper block continues vertically to
    here.
    _ (underline char): means that the previous structure continues
    here.
    - (dash): empty space

    Jul 08 2022
    Saeed
    """


class Structure:
    """read the struct file"""
    def __init__(self) -> None:
        self.strcut: str = sys.argv[1]

    def mk_matrix(self) -> None:
        """make a matrix out of the blocks symbols"""

    def read_struct(self) -> tuple[dict, dict]:
        """read the strut file"""
        f: typing.IO  # a string to save file
        line: str  # a string to save lines of the strcut file
        bed_count: int = 0  # to count lines in the matrix of bolcks
        symbole_dict: dict[str, str] = {}  # dict to save name and symb
        block_dict: dict[int, list[str]] = {}  # dict to save matrix

        with open(self.strcut, 'r') as f:
            while True:
                line = f.readline()
                if line.strip().startswith('#'):
                    pass
                elif line.strip().startswith('!'):
                    sym, fname = self.get_files(line.strip())
                    symbole_dict[sym] = fname
                elif line.strip():
                    m_list = self.get_matrix(line.strip())
                    block_dict[bed_count] = m_list
                    bed_count += 1
                if not line:
                    break
        self.check_dicts(symbole_dict, block_dict)
        return symbole_dict, block_dict

    def get_files(self, line: str) -> tuple[str, str]:
        """check the files name and if they are not empty"""
        # Drop ! from beginning
        line = re.sub('!', '', line)
        # Remove whit spaces
        line = re.sub(r'\s+', '', line)
        sym, fname = line.split("=")
        self.check_files(fname)
        return sym, fname

    def check_files(self, fname: str) -> None:
        """check if the fname exist and not empty"""
        if not os.path.isfile(fname):
            exit(f'ERROR: "{fname}" does not exist!!\n')
        if not os.path.getsize(fname) > 0:
            exit(f'ERROR: "{fname}" is empty!!\n')

    def get_matrix(self, line: str) -> list[str]:
        """read the matrix section of the struct file"""
        _sym_mat: list[str]  # A list to return sequence in line
        if ' ' in line:
            print(f'WAR: whitespace in the "{self.strcut}" in line: '
                  f'"{line}", it removed!\n')
            line = re.sub(r'\s+', '', line)
        _sym_mat = [item for item in line]
        return _sym_mat

    def check_dicts(self,
                    sym: dict[str, str],
                    block: dict[int, list[str]]) -> tuple[dict, dict]:
        """check all the symbols have a file defeind with them"""
        e_flag: bool = False  # To check all the typo in the input file
        for _, row in block.items():
            for i in range(len(row)):
                if row[i].isalpha():
                    if not row[i] in sym.keys():
                        print(f'ERROR: symbole "{row[i]}" is not defined\n')
                        e_flag = True
                elif row[i] not in ['-', '_', '|']:
                    print(f'ERROR: symbole "{row[i]}" is not defined\n')
                    e_flag = True
        if e_flag:
            exit(f'Mistakes in the "{self.strcut}"')


super_str = Structure()
super_str.read_struct()
