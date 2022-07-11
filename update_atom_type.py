import typing
import pandas as pd
import read_lmp_data as mlmp


class Doc:
    """update atom type for combination of data
    input a dict in the following form:
    {'A': 'decane.data', 'B: 'water.data'}
    """


class UpdateType:
    """updata types"""
    def __init__(self, files: dict[str, str]) -> None:
        self.files = files
        del files
        self.up_files = self.update_atom_type()

    def update_atom_type(self) -> dict[str, dict[str, typing.Any]]:
        atom_indent: int = 0  # to increase the type of each file
        bond_indent: int = 0  # to increase the type of each file
        angle_indent: int = 0  # to increase the type of each file
        dihedral_indent: int = 0  # to increase the type of each file
        up_dict: dict[str, dict[str, typing.Any]] = dict()
        for i, (k, v) in enumerate(self.files.items()):
            up_dict[k] = dict()
            up_dict[k]['fname'] = v
            read_data = mlmp.ReadData(v)
            print(f'{self.__class__.__name__}:\n'
                  f'\tUpdating: {v}\n')
            if i == 0:
                atom_indent += read_data.NAtomTyp
                bond_indent += read_data.NBondTyp
                angle_indent += read_data.NAngleTyp
                dihedral_indent += read_data.NDihedralTyp
            else:
                if read_data.NAtomTyp > 0:
                    read_data.Atoms_df = self._update_type(
                        read_data.Atoms_df, atom_indent)
                    atom_indent += read_data.NAtomTyp
                else:
                    exit(f'{self.__class__.__name__}:\n'
                         f'\tERROR:  ZERO atom type in "{v}"')
                if read_data.NBondTyp > 0:
                    read_data.Bonds_df = self._update_type(
                        read_data.Bonds_df, bond_indent)
                    bond_indent += read_data.NBondTyp
                else:
                    pass
                if read_data.NAngleTyp > 0:
                    read_data.Angles_df = self._update_type(
                        read_data.Angles_df, angle_indent)
                    angle_indent += read_data.NAngleTyp
                else:
                    pass
                if read_data.NDihedrals > 0:
                    read_data.Dihedrals_df = self._update_type(
                        read_data.Dihedrals_df, dihedral_indent)
                    dihedral_indent += read_data.NDihedralTyp
                else:
                    pass
            up_dict[k]['data'] = read_data
        return up_dict

    def _update_type(self, df: pd.DataFrame, indent: int) -> pd.DataFrame:
        df['typ'] += indent
        return df
