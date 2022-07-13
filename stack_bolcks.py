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
        row_list = self.stack_atoms()
        self.stack_rows(row_list)
        del block, bs, row_list

    def stack_atoms(self) -> list[pd.DataFrame]:
        """stack atoms with updating id of atoms and also mol number"""
        # Declear variables
        col: int      # Column of the block data in struct file
        row: int      # Row of the block data in struct file
        v: list(str)  # List of symbols of the data file
        item: str     # An item of v

        self.VACUME: int = 2  # Space between blocks

        x_indent: float = 0.0  # Shift in x column based on previous block

        max_x: float = 0  # max in x column after updating DataFrame
        max_id: int = 0   # max of id column after updatinf DataFrame
        max_mol: int = 0  # max of mol column after updatinf DataFrame

        _df: pd.DataFrame  # Temperary dataframe
        raw_df: pd.DataFrame  # A Dataframe to keep all the row made DF
        row_list: list[pd.DataFrame] = []  # A list to keep the raws DF
        for row, (_, v) in enumerate(self.block.items()):
            df_list: list[pd.DataFrame] = []  # To append DataFrames
            for col, item in enumerate(v):
                _df = self.bs.system[item]['data'].Atoms_df.copy()
                if col == 0:
                    df_list.append(_df)
                    max_x = np.max(_df['x'])
                    max_id = np.max(_df['atom_id'])
                    max_mol = np.max(_df['mol'])
                    del _df
                else:
                    x_indent = max_x + self.VACUME
                    _df['x'] += x_indent
                    _df['atom_id'] += max_id
                    _df['mol'] += max_mol
                    max_x = np.max(_df['x'])
                    max_id = np.max(_df['atom_id'])
                    max_mol = np.max(_df['mol'])
                    df_list.append(_df)
                    del _df
            raw_df = pd.concat(df_list, ignore_index=True,  axis=0)
            row_list.append(raw_df)
            del raw_df, df_list
        return row_list

    def stack_rows(self, row_list: list[pd.DataFrame]) -> None:
        max_z: float = 0  # max in z column after updating DataFrame
        z_indent: float = 0  # Shift in z column based on previous block
        x_ave_center: float  # x of center_of_mass of bottom layer
        y_ave_center: float  # y of center_of_mass of bottom layer
        x_ave: float  # x center_of_mass of each layer
        y_ave: float  # y center_of_mass of each layer
        Atoms_df: pd.DataFrame  # To return data
        for i, item in enumerate(row_list):
            if i == 0:
                max_id = np.max(item['atom_id'])
                max_mol = np.max(item['mol'])
                max_z = np.max(item['z'])
                x_ave_center = 0.5*(np.max(item['x'] - np.min(item['x'])))
                y_ave_center = 0.5*(np.max(item['y'] - np.min(item['y'])))
            else:
                z_indent = max_z + self.VACUME
                item['z'] += z_indent
                item['mol'] += max_mol
                item['atom_id'] += max_id
                max_id = np.max(item['atom_id'])
                max_mol = np.max(item['mol'])
                max_z = np.max(item['z'])
                x_ave = 0.5*(np.max(item['x'] - np.min(item['x'])))
                y_ave = 0.5*(np.max(item['y'] - np.min(item['y'])))
                item['x'] -= (x_ave-x_ave_center)
                item['y'] -= (y_ave-y_ave_center)
        Atoms_df = pd.concat(row_list, ignore_index=True,  axis=0)
        self.Atoms_df = Atoms_df
