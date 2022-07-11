import numpy as np
import pandas as pd


class Doc:
    """stack atom in Z loop
    first in x and then in y direction"""


class stack_data:
    """put atoms side to side"""
    def __init__(self,
                 block: pd.DataFrame,  # Block information
                 bs: dict[str, dict[str, str]]  # Inforamtion for all systems
                 ) -> None:
        self.block = block
        self.bs = bs
        del block, bs
        # print(self.bs.system['W']['data'].Atoms_df)
        self.stack_atoms()

    def stack_atoms(self):
        """stack atoms with updating id of atoms and also mol number"""
        row: int  # Row of the block data in struct file
        v: list(str)  # List of symbols of the data file
        col: int  # Column of the block data in struct file
        item: str  # An item from v

        for row, (_, v) in enumerate(self.block.items()):
            for col, item in enumerate(v):
                print(row, col, item)
                # print(self.bs.system[item]['data'].Atoms_df)
