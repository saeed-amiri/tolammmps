import sys
import os
import re
import typing
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

    def read_struct(self) -> None:
        """read the strut file"""
        bed_count: int = 0  # to count lines in the matrix of bolcks
        f: typing.IO  # a string to save file
        line: str  # a string to save lines of the strcut file
        symbole_dict: dict[str, str] = {} # dict to save name and symb
        mat_dict: dict[str, list[str]] = {} # dict to save matrix

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
                    mat_dict[bed_count] = m_list
                    bed_count += 1
                if not line:
                    break
        print(mat_dict)

    def get_files(self, line: str) -> tuple[str, str]:
        """check the files name and if they are not empty"""
        # Drop ! from beginning
        line = re.sub('!', '', line)
        # Remove whit spaces
        line = re.sub(r'\s+', '', line)
        sym, fname = line.split("=")
        return sym, fname
    
    def get_matrix(self, line: str) -> None:
        """read the matrix section of the struct file"""
        _sym_mat: list[str]  # A list to return sequence in line
        if ' ' in line: 
            print(f'WAR: whitespace in the "{self.strcut}" in line: '
                  f'"{line}", it removed!\n')
            line = re.sub(r'\s+', '', line)
        _sym_mat = [item for item in line]
        return _sym_mat


super_str = Structure()
super_str.read_struct()
