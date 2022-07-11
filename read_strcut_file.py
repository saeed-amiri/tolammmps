import os
import re
import sys
import typing

class Structure:
    """read the struct file"""
    def __init__(self) -> None:
        self.fname: str = sys.argv[1]

    def mk_block(self) -> None:
        """make a matrix out of the blocks symbols"""
        self.check_files(self.fname)
        print(f'{self.__class__.__name__}:\n'
              f'\treading "{self.fname}"')
        self.symbols, self.block = self.read_struct()

    def read_struct(self) -> tuple[dict, dict]:
        """read the strut file"""
        f: typing.IO  # a string to save file
        line: str  # a string to save lines of the strcut file
        bed_count: int = 0  # to count lines in the matrix of bolcks
        symbole_dict: dict[str, str] = {}  # dict to save name and symb
        block_dict: dict[int, list[str]] = {}  # dict to save matrix

        with open(self.fname, 'r') as f:
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
            exit(f'{self.__class__.__name__}:\n'
                 f'\tERROR: "{fname}" does not exist!!\n')
        if not os.path.getsize(fname) > 0:
            exit(f'{self.__class__.__name__}:\n'
                 f'\tERROR: "{fname}" is empty!!\n')

    def get_matrix(self, line: str) -> list[str]:
        """read the matrix section of the struct file"""
        _sym_mat: list[str]  # A list to return sequence in line
        if ' ' in line:
            print(f'{self.__class__.__name__}:\n'
                  f'\t"{self.fname}" -> WARRNING: whitespace in the in line: '
                  f'"{line}", it is removed!\n')
            line = re.sub(r'\s+', '', line)
        _sym_mat = [item for item in line]
        return _sym_mat

    def check_dicts(self,
                    sym: dict[str, str],
                    block: dict[int, list[str]]) -> None:
        """check if all the symbols have a file defeind with them"""
        e_flag: bool = False  # To check all the typo in the input file
        for _, row in block.items():
            for i in range(len(row)):
                if row[i].isalpha():
                    if not row[i] in sym.keys():
                        print(f'{self.__class__.__name__}:\n'
                              f'\tERROR: "{self.fname}" -> symbole "{row[i]}"'
                              f' is not defined\n')
                        e_flag = True
                elif row[i] not in ['-', '_', '|']:
                    print(f'{self.__class__.__name__}:\n'
                          f'\tERROR: "{self.fname}" -> symbole "{row[i]}" is '
                          f'not defined\n')
                    e_flag = True
        if e_flag:
            exit(f'Mistake(s) in the "{self.fname}"')


if __name__ == "__main__":
    super_str = Structure()
    super_str.mk_block()
    print(super_str.block)
