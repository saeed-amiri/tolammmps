import os
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
        self.rotate_ax(df)

    def rotate_ax(self, df: pd.DataFrame) -> None:
        """rotate alnog x axis"""
        Nmols = np.max(df['mol'])
        self.random_tetha(Nmols)

    def random_tetha(self, Nmols: int) -> list[float]:
        """return a list of random angles, with length Nmols"""
        theta_list: list[float] = np.random.normal(size=Nmols)
        theta_list = [item*np.pi/2 for item in theta_list]
        print(theta_list)

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


if __name__ == '__main__':
    fname = sys.argv[1]
    data = rlmp.ReadData(fname)
    rot = Rotate(data.Atoms_df)
