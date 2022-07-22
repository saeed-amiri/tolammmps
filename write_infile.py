import re
import sys
import typing
import datetime
import pandas as pd
from datetime import date


class Doc:
    """write LAMMPS input file based on the parameter"""


class WriteIn:
    """write input file"""
    def __init__(self) -> None:
        self.fname = 'in.lmp'  # Output file
        self.write_in()

    def write_in(self) -> None:
        """call all the method"""
        with open(self.fname, 'w') as f:
            self.write_definitions(f)
            self.write_computes(f)
    
    def write_definitions(self, f) -> None:
            self.write_header(f)
            self.write_info(f)
            self.write_init(f)
            self.write_gpu(f)
            self.write_include(f)
            self.write_minimization(f)
    
    def write_computes(self, f) -> None:
        """write computation section"""
        f.write(f"#{'Computation':.^85}\n")
        f.write(f'\n')
        f.write(f"#{'':_^85}\n")

    def write_minimization(self, f) -> None:
        """write minimiation section"""
        f.write(f"#{'Minimization':.^85}\n")
        f.write(f'\tminimize 1.0e-3 1.0e-3 1000 10000\n')
        f.write(f'\twrite_data min.data\n')
        f.write(f'\treset_timestep 0\n')
        f.write(f"#{'':_^85}\n")
        f.write(f'\n')


    def write_include(self, f: typing.TextIO) -> None:
        """write commands to include or load files"""
        f.write(f"#{'Include':.^85}\n")
        f.write(f'\tread_data \n')
        f.write(f'\tinclude \n')
        f.write(f"#{'':_^85}\n")
        f.write(f'\n')


    def write_init(self, f: typing.TextIO) -> None:
        """write initiation of the simulation"""
        f.write(f"#{'Initiation':.^85}\n")
        f.write(f'\tunits \n')
        f.write(f'\ttimestep \n')
        f.write(f'\tatom_style \n')
        f.write(f'\tneigh_modify every 1 delay 0 check yes\n')
        f.write(f'\tspecial_bonds lj/coul 0.0 0.0 0.5\n')
        f.write(f'\tkspace_style pppm 0.0001\n')
        f.write(f"#{'':_^85}\n")
        f.write(f'\n')


    def write_gpu(self, f: typing.TextIO) -> None:
        f.write(f"#{'GPU specification':.^85}\n")
        f.write(f'\tpackage gpu 1\n')
        f.write(f'\tpackage gpu 1 neigh yes\n')
        f.write(f'\tpackage gpu 1 binsize 1\n')
        f.write(f"#{'':_^85}\n")
        f.write(f'\n')
    def write_info(self, f: typing.TextIO) -> None:
        """write info about input file"""
        fname = f'fname: {sys.argv[1]}'
        f.write(f'#{"Info:":^85}\n')
        f.write(f'#{fname:^85}\n')
        f.write(f"#{'':_^85}\n")
        f.write(f'\n')


    def write_header(self, f: typing.TextIO) -> None:
        """write header of the file"""
        v: str = 'LAMMPS'
        today: datetime.date = date.today()
        code: str = f'class: {self.__class__.__name__}'
        i_file: str = f'script: {__file__}'
        f.write(f"#{'':,^85}\n")
        f.write(f"#{'':`^85}\n")
        f.write(f'#{v:^85}\n')
        f.write(f'#\n')
        f.write(f"#{i_file:^85}\n")
        f.write(f"#{code:^85}\n")
        f.write(f'#{str(today):^85}\n')
        f.write(f"#{'':`^85}\n")


if __name__ == "__main__":
    w = WriteIn()
