import os
from pprint import pprint
import sys
import pandas as pd
import numpy as np
import read_lmp_data as rlmp
import write_lmp as wlmp


class Doc:
    """rotate molecules
    Rotate molecules around thier axis randomly
    Input:
        LAMMPS data file
    Output:
        LAMMPS data file
    """


class Rotate:
    """get the data file and rotate along axis"""
    def __init__(self, df: pd.DataFrame) -> None:
        self.rotate(df)

    def rotate(self, df: pd.DataFrame) -> None:
        """call functions"""
        # self.df = self.same_rotate_ax(df)
        self.df = self.rotate_com(df)
        self.bent_mol(df)

    def rotate_com(self, df: pd.DataFrame) -> pd.DataFrame:
        """rotate each molecule around its center of mass (aom)"""
        mol_list = df['mol'].copy()
        mol_set = set(mol_list)
        Nmols: int = len(mol_set)
        # Make a list of random angles to rotate molecules:
        theta_list: list[float]  # Random angles for each molecule
        theta_list = self.random_tetha(Nmols)
        com_list: list[tuple[float, float, float]]  # COM of each molecule
        com_list = self.get_com(df, mol_set)
        x_list = []
        y_list = []
        z_list = []
        for i, atom in df.iterrows():
            mol_index = atom['mol'] - 1
            x_ = atom['x'] - com_list[mol_index][0]
            y_ = atom['y'] - com_list[mol_index][1]
            z_ = atom['z'] - com_list[mol_index][2]
            p_x, p_y, p_z = self.y_rotation(theta_list[mol_index], x_, y_, z_)
            x_list.append(p_x + com_list[mol_index][0])
            y_list.append(p_y + com_list[mol_index][1])
            z_list.append(p_z + com_list[mol_index][2])
        df['x'] = x_list
        df['y'] = y_list
        df['z'] = z_list
        return df

    def bent_mol(self, df: pd.DataFrame) -> pd.DataFrame:
        """benting the moles by transformin them to spherical coords"""
        mol_list = df['mol'].copy()
        mol_set = set(mol_list)
        com_list: list[tuple[float, float, float]]  # COM of each molecule
        com_list = self.get_com(df, mol_set)
        x_list = []
        y_list = []
        z_list = []
        for i, atom in df.iterrows():
            mol_index = atom['mol'] - 1
            x_ = atom['x'] - com_list[mol_index][0]
            y_ = atom['y'] - com_list[mol_index][1]
            z_ = atom['z'] - com_list[mol_index][2]
            p_x, p_y, p_z = self.cartisain_to_spherical(x_, y_, z_)
            x_list.append(p_x)
            y_list.append(p_y)
            z_list.append(p_z)
        df['x'] = x_list
        df['y'] = y_list
        df['z'] = z_list
        return df

    def get_com(self,
                df: pd.DataFrame,
                mol_set: set[float]) -> list[tuple[float, float, float]]:
        """get the center of mass of each molecule"""
        # A list to append the COM of each molecule:
        com_list: list[tuple[float, float, float]] = []
        df_i: pd.DataFrame  # temp DataFrame
        for mol in mol_set:
            df_i = df.loc[df['mol'] == mol].copy()
            com_list.append(self.return_com(df_i))
            del df_i
        return com_list

    def return_com(self, df: pd.DataFrame) -> tuple[float, float, float]:
        """return mean of each columns *x or y or z"""
        x_com: float = df['x'].mean()
        y_com: float = df['y'].mean()
        z_com: float = df['z'].mean()
        return x_com, y_com, z_com

    def same_rotate_ax(self, df: pd.DataFrame) -> pd.DataFrame:
        """rotate alnog x axis"""
        Nmols: int = np.max(df['mol'])
        theta_list: list[float] = self.random_tetha(Nmols)
        rot_x: list[float] = []  # To append the new x
        rot_y: list[float] = []  # To append the new y
        rot_z: list[float] = []  # To append the new z
        for x, y, z, mol in zip(df['x'], df['y'], df['z'], df['mol']):
            x_, y_, z_ = self.x_rotation(theta_list[mol-1], x, y, z)
            rot_x.append(x_)
            rot_y.append(y_)
            rot_z.append(z_)
        df['x'] = rot_x
        df['y'] = rot_y
        df['z'] = rot_z
        del rot_x, rot_y, rot_z
        return df

    def random_tetha(self, Nmols: int) -> list[float]:
        """return a list of random angles, with length Nmols"""
        theta_list: list[float] = np.random.normal(size=Nmols)
        theta_list = [item*np.pi/6 for item in theta_list]
        return theta_list

    def x_rotation(self,
                   theta: float,
                   x: float,
                   y: float,
                   z: float) -> tuple[float, float, float]:
        """calculate the rotation matrix along x axis"""
        x_: float = x
        y_: float = y*np.cos(theta)-z*np.sin(theta)
        z_: float = y*np.sin(theta)+z*np.cos(theta)
        return x_, y_, z_

    def y_rotation(self,
                   theta: float,
                   x: float,
                   y: float,
                   z: float) -> tuple[float, float, float]:
        """calculate the rotation matrix along y axis"""
        x_: float = x*np.cos(theta)+z*np.sin(theta)
        y_: float = y
        z_: float = z*np.cos(theta)-x*np.sin(theta)
        return x_, y_, z_

    def z_rotation(self,
                   theta: float,
                   x: float,
                   y: float,
                   z: float) -> tuple[float, float, float]:
        """calculate the rotation matrix along z axis"""
        x_: float = x*np.cos(theta)-y*np.sin(theta)
        y_: float = y*np.sin(theta)+y*np.cos(theta)
        z_: float = z
        return x_, y_, z_

    def cartisain_to_spherical(self, x, y, z) -> tuple[float, float, float]:
        """Conversion between Cylindrical and Cartesian Coordinates"""
        r: float = np.sqrt(x*x+y*y+z*z)
        theta: float = np.arctan2(y, x)
        phi: float = np.arctan2(np.sqrt(x*x+y*y), z)
        return r, theta, phi


if __name__ == '__main__':
    fname = sys.argv[1]
    data = rlmp.ReadData(fname)
    rot = Rotate(data.Atoms_df)
    data.Atoms_df = rot.df
    wrt = wlmp.WriteLmp(data)
    wrt.write_lmp()
