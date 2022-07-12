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
        item: str     # An item from v

        x_indent: float = 0.0  # Shift in x column based on previous block
        y_indent: float = 0.0  # Shift in y column based on previous block
        z_indent: float = 0.0  # Shift in z column based on previous block

        x_vacume: int = 2  # Space between blocks

        max_x: float  # max in x column after updating DataFrame
        max_y: float  # max in y column after updating DataFrame
        max_z: float  # max in z column after updating DataFrame
        max_id: int   # max of id column after updatinf DataFrame
        max_mol: int  # max of mol column after updatinf DataFrame

        _df: pd.DataFrame  # Temperary dataframe
        Atoms_df: pd.DataFrame  # To return data

        for row, (_, v) in enumerate(self.block.items()):
            for col, item in enumerate(v):
                if row == 0:
                    if col == 0:
                        Atoms_df = self.bs.system[item]['data'].Atoms_df
                        max_x = np.max(Atoms_df['x'])
                        max_id = np.max(Atoms_df['atom_id'])
                        max_mol = np.max(Atoms_df['mol'])
                    else:
                        _df = self.bs.system[item]['data'].Atoms_df
                        _df['x'] += max_x
                        _df['x'] += x_vacume
                        _df['atom_id'] += max_id
                        _df['mol'] += max_mol
                        Atoms_df = pd.concat([Atoms_df, _df])
                        max_x = np.max(Atoms_df['x'])
                        max_id = np.max(Atoms_df['atom_id'])
                        max_mol = np.max(Atoms_df['mol'])
                        del _df
                else:
                    max_y = np.max(Atoms_df['y'])
                    if col == 0:
                        _df = self.bs.system[item]['data'].Atoms_df
                        _df['x'] += max_x
                        _df['x'] += x_vacume
                        _df['atom_id'] += max_id
                        _df['mol'] += max_mol
                        Atoms_df = pd.concat([Atoms_df, _df])
                        del _df
