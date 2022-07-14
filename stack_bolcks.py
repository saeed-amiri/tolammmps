import numpy as np
import pandas as pd


class Doc:
    """stack atom in Z loop
    first in x and then in y direction"""


class UpdateAtoms:
    """put atoms side to side"""
    def __init__(self,
                 block: pd.DataFrame,  # Block information
                 bs: dict[str, dict[str, str]]  # Inforamtion for all systems
                 ) -> None:

        self.block = block
        self.bs = bs
        self.VACUME: int = 2  # Space between blocks
        self.stack_atoms()
        del block, bs

    def stack_atoms(self) -> None:
        """stack atoms along axises"""
        print(f'{self.__class__.__name__}:\n'
              f'\tUpdating: Atoms\n')
        row_list = self.stack_x()
        self.stack_rows(row_list)

    def stack_x(self) -> list[pd.DataFrame]:
        """stack atoms with updating id of atoms and also mol number
        this function stack data along x-axis"""
        # Declear variables
        col: int      # Column of the block data in struct file
        row: int      # Row of the block data in struct file
        v: list[str]  # List of symbols of the data file
        item: str     # An item of v
        raw_df: pd.DataFrame  # A Dataframe to keep all the row made DF
        row_list: list[pd.DataFrame] = []  # A list to keep the raws DF
        _df: pd.DataFrame  # Temperary dataframe

        x_indent: float = 0.0  # Shift in x column based on previous block
        max_x: float = 0  # max in x column after updating DataFrame
        max_id: int = 0   # max of id column after updatinf DataFrame
        max_mol: int = 0  # max of mol column after updatinf DataFrame
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
        """stack atoms with updating id of atoms and also mol number
        this function stack data from self.stack_atoms along y-axis"""
        # Declear variables
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
        self.Atoms_df = pd.concat(row_list, ignore_index=True,  axis=0)


class UpdateBond:
    """updating bonds dataframe in and stack them in one DataFrame"""
    def __init__(self,
                 block: pd.DataFrame,  # Block information
                 bs: dict[str, dict[str, str]]  # Inforamtion for all systems
                 ) -> None:
        self.block = block
        self.bs = bs
        self.update_bonds()
        del bs, block

    def update_bonds(self):
        _df: pd.DataFrame  # A temprarry df to store them
        df_list: list[pd.DataFrame] = []  # To append all the DFs
        prev_natoms: int = 0  # Number of atoms up to current
        prev_nbonds: int = 0  # Number of bonds up to current
        row: int  # counting the rows, not used here just for clairity
        col: int  # counting the cols, not used here just for clairity
        for row, (_, v) in enumerate(self.block.items()):
            for col, item in enumerate(v):
                _df = self.bs.system[item]['data'].Bonds_df.copy()
                prev_nbonds += self.bs.system[item]['data'].NBonds
                _df['ai'] += prev_natoms
                _df['aj'] += prev_natoms
                df_list.append(_df)
                prev_natoms += self.bs.system[item]['data'].NAtoms
        self.Bonds_df = pd.concat(df_list, ignore_index=True,  axis=0)
        self.Bonds_df.index += 1
        if prev_nbonds != len(self.Bonds_df):
            exit(f'\tERROR!: Problem in number of bonds\n'
                 f'\tnumber of total bonds: {prev_nbonds}'
                 f' != {len(self.Bonds_df)} number of calculated bonds')

class StackData(UpdateAtoms, UpdateBond):
    """stack all the DataFrame together"""
    def __init__(self,
                 block: pd.DataFrame,
                 bs: dict[str, dict[str, str]]
                 ) -> None:
        UpdateAtoms.__init__(self, block, bs)
        UpdateBond.__init__(self, block, bs)
