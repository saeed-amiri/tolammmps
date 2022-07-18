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
        param_dict: dict[str, list[str]] = dict()  # to save lines from fname
        data_dict: dict[int, list[str]] = dict()  # To save lines for each file
        line: str  # string to save read lines
        append_flag: bool = False  # Flage to see the data files
        file_count: int = -1  # file_id to lable each block
        with open(self.fname, 'r') as f:
            while True:
                line = f.readline()
                if line.strip().startswith('{'):
                    file_count += 1
                    append_flag = True
                    data_dict[file_count] = []
                elif line.startswith('}'):
                    append_flag = False
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
                symb = self.return_item_value(item)
            elif item.startswith('file'):
                fname = self.return_item_value(item)
            elif item.startswith('style'):
                style = self.return_item_value(item)

    def return_item_value(self, item: str) -> str:
        """split and strip the vale of each line"""
        return item.split('=')[1].strip()


if __name__ == '__main__':
    param = ReadParameter()
