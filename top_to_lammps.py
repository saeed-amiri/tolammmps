import os, sys, re
from pprint import pprint
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
        self.FLAG_list = []
        self.FORMAT_list = []
        pass

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
        free_dict = [dict() for item in self.FORMAT_list]
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

class READTOP (TOP):
    """%FLAG POINTERS                                                                  
    %FORMAT(10I8)                                                                   
    6702       7    5350       0       0       0       0       0       0       0
    8468    3062       0       0       0       3       0       0      10       1
       0       0       0       0       0       0       0       1       8       0
       0

    The cards are space-parsed in some cases, i.e. names of atoms, mols, residues
    """
    def __init__(self) -> None:
        super().__init__()
        self.read_file()

    def get_data(self)->None:
        self.mk_modify()

    def mk_modify(self)-> None:
        """modify the FLAG dict"""
        self.set_flage()

    def set_flage(self) -> list:
        """set True and False flag to cards' name"""
        for key in self.FLAG: self.FLAG[key]['flag']=False

        pprint(self.FLAG)
if __name__=="__main__":
    TOPFILE = "test3.top"
    top = READTOP()
    top.get_data()
    # top.read_file()

    
