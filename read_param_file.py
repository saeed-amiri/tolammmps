import re
import os
import sys
import pandas as pd
from pprint import pprint


class Doc:
    """read the 'param.data': parameter file to write LAMMPS
    parameter.lmp file
        Input:
            param.data, which should contain:
                pair style
                style coefficients
                bond, angle, dihedral styles, and parameters
        Out:
            pd.Dataframe contains all the parameters of interactions.
    """


class ReadParameter:
    """read parmater file"""
    def __init__(self) -> None:
        self.fname: str = 'param.data'
        self.read_param()

    def read_param(self) -> None:
        """read parameter file"""
        data_dict: dict[int, list[str]] = self.get_data_block()
        self.process_data(data_dict)

    def get_data_block(self) -> dict[int, list[str]]:
        line: str  # String to save read lines
        file_count: int = -1  # "file_id" to lable each block
        append_flag: bool = False  # Flage to see the data files
        data_dict: dict[int, list[str]] = dict()  # To save lines for each file

        with open(self.fname, 'r') as f:
            while True:
                line = f.readline()
                if line.strip().startswith('{'):
                    file_count += 1
                    append_flag = True
                    data_dict[file_count] = []
                elif line.startswith('}'):
                    append_flag = False
                elif line.startswith('#'):
                    pass
                elif append_flag:
                    data_dict[file_count].append(line.strip())
                if not line:
                    break
        return data_dict

    def process_data(self, data: dict[int, list[str]]) -> None:
        """process each block of data"""
        for file_id, params in data.items():
            self.set_block_param(params)

    def set_block_param(self, parameters: list[str]) -> None:
        """set parameters for each file"""
        for item in parameters:
            if item.startswith('symb'):
                symb = self.return_item_value(item, 1)
            elif item.startswith('file'):
                fname = self.return_item_value(item, 1)
            elif item.startswith('style'):
                style = self.return_item_value(item, 1)
            elif item.startswith('@'):
                self.return_atom(item)

    def return_atom(self, item: str) -> list[str]:
        """return a list contains info of each atom"""
        atom_info = item.strip("@")
        atom_type = self.return_item_value(atom_info, 0)
        atom_param = self.return_item_value(atom_info, 1)
        print(self.set_atom_atrr(atom_type, atom_param))

    def set_atom_atrr(self, atom_type: int, atom: str) -> dict[str, list[str]]:
        """set the properties of the atom, as a attributes to thier 
        name"""
        # drop brackets
        atom = re.sub(r"[\([{})\]]", "", atom)
        atom_attr = atom.split(',')
        atom_dict: dict[int, str] = dict()  # To save info for each atom
        for item in atom_attr:
            attrs = self.return_item_value(item, 0)
            value = self.return_item_value(item, 1)
            atom_dict[attrs] = value
        return atom_dict

    def return_item_value(self, item: str, loc: int) -> str:
        """split and strip the vale of each line"""
        return item.split('=', 1)[loc].strip()



if __name__ == '__main__':
    param = ReadParameter()
