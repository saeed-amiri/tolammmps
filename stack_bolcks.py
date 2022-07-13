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
        self.stack_atoms()
        # print(self.bs.system['W']['data'].Atoms_df)

    def stack_atoms(self):
        """stack atoms with updating id of atoms and also mol number"""
        # Declear variables
        col: int      # Column of the block data in struct file
        row: int      # Row of the block data in struct file
        v: list(str)  # List of symbols of the data file
        item: str     # An item of v

        x_vacume: int = 2  # Space between blocks

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
        df_list: list[pd.DataFrame] = [] # To append DataFrames
        df_dict = dict()
        for row, (_, v) in enumerate(self.block.items()):
            for col, item in enumerate(v):
                _df = self.bs.system[item]['data'].Atoms_df.copy()
                if row == 0:
                    max_y = np.max(_df['y'])
                else:
                    _df['y'] += max_y
                    _df['y'] += x_vacume
                if col == 0:
                    df_list.append(_df)
                    max_x = np.max(_df['x'])
                    max_id = np.max(_df['atom_id'])
                    max_mol = np.max(_df['mol'])
                    del _df
                else:
                    _df['x'] += max_x
                    _df['x'] += x_vacume
                    _df['atom_id'] += max_id
                    _df['mol'] += max_mol
                    max_x = np.max(_df['x'])
                    max_id = np.max(_df['atom_id'])
                    max_mol = np.max(_df['mol'])
                    df_list.append(_df)
                    del _df

        Atoms_df = pd.concat(df_list, ignore_index=True,  axis=0)
        # print(Atoms_df)
        self.Atoms_df = Atoms_df
