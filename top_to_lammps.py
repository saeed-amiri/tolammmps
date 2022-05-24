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
        self.read_card()

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
                            for key in self.FLAG.keys():self.FLAG[key]['flag'] = False
                            flag = line.split('FLAG')[1].strip()
                            if flag in self.FLAG.keys(): 
                                self.FLAG[flag]['flag'] = True
                        elif line.startswith("FORMAT"): pass
                else:
                    # reading data for each card
                    for key in self.FLAG.keys():
                        if self.FLAG[key]['flag']:
                            self.FLAG[key]['data'].append(line)
                if not line: break

if __name__== "__main__":
    TOPFILE = "test3.top"
    top = READTOP()
    top.get_data()
    # top.read_file()
    # for key in top.FLAG.keys(): print(key)
    pprint(top.FLAG['ATOM_NAME'])

    
