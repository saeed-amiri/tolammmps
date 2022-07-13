import numpy as np
import pandas as pd
from sqlalchemy import column

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
        self.stack_atoms()
        # print(self.bs.system['W']['data'].Atoms_df)

    def stack_atoms(self):
        """stack atoms with updating id of atoms and also mol number"""
        # Declear variables
        col: int      # Column of the block data in struct file
        row: int      # Row of the block data in struct file
        v: list(str)  # List of symbols of the data file
        item: str     # An item of v

        VACUME: int = 2  # Space between blocks

        x_indent: float = 0.0  # Shift in x column based on previous block
        y_indent: float = 0.0  # Shift in y column based on previous block
        z_indent: float = 0.0  # Shift in z column based on previous block

        max_x: float = 0 # max in x column after updating DataFrame
        max_y: float = 0 # max in y column after updating DataFrame
        max_z: float = 0 # max in z column after updating DataFrame
        max_id: int = 0   # max of id column after updatinf DataFrame
        max_mol: int = 0  # max of mol column after updatinf DataFrame

        _df: pd.DataFrame  # Temperary dataframe
        Atoms_df: pd.DataFrame  # To return data
        raw_df: pd.DataFrame  # A Dataframe to keep all the row made DF
        raw_list: list[pd.DataFrame] = []  # A list to keep the raws DF
        for row, (_, v) in enumerate(self.block.items()):
            df_list: list[pd.DataFrame] = [] # To append DataFrames
            for col, item in enumerate(v):
                _df = self.bs.system[item]['data'].Atoms_df.copy()
                if col == 0:
                    df_list.append(_df)
                    max_x = np.max(_df['x'])
                    max_id = np.max(_df['atom_id'])
                    max_mol = np.max(_df['mol'])
                    del _df
                else:
                    x_indent = max_x + VACUME
                    _df['x'] += x_indent
                    _df['atom_id'] += max_id
                    _df['mol'] += max_mol
                    max_x = np.max(_df['x'])
                    max_id = np.max(_df['atom_id'])
                    max_mol = np.max(_df['mol'])
                    df_list.append(_df)
                    del _df
            raw_df = pd.concat(df_list, ignore_index=True,  axis=0)
            raw_list.append(raw_df)
            del raw_df, df_list
        for i, item in enumerate(raw_list):
            if i == 0:
                max_id = np.max(item['atom_id'])
                max_mol = np.max(item['mol'])
                max_y = np.max(item['y'])
            else:
                item['atom_id'] += max_id
                item['mol'] += max_mol
                y_indent = max_y + VACUME
                item['y'] += y_indent
                max_id = np.max(item['atom_id'])
                max_mol = np.max(item['mol'])
                max_y = np.max(item['y'])
        Atoms_df = pd.concat(raw_list, ignore_index=True,  axis=0)
        self.Atoms_df = Atoms_df
