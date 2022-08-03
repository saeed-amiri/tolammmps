import read_strcut_file as struct
import update_atom_type as mupdate
import stack_bolcks as mstack
import write_lmp as lmpwrt
import write_param as parwrt


class Doc:
    """This code combines data files and returns a superstructure in
    LAMMPS full atom style.   The input file is <system>.struct file
    which contains the names and path of the  files and also a matrix
    that shows the order of the superstructure:

   ex. of an input file:

    First, it should have a symbol for each file and its path.
    The symbol is better to be an upper letter (for now, it does not
    recognize lower/upper case of the letter):

        ! D=decane.data
        ! W=water.data
        ! S=sio2.data
    Also the second axis of stacking can be choosen between y and z
    the defaulf value is z:
    axis = y
    Then it should have a matrix that shows how you want to build your
    superstructure:
        DWD
        _S_
        DWD
    which stack the blocks as: "decane water decane" most top level,
    SiO2 will stack in the second layer and "decane water decane" in
    the lowest layer. The atoms id will be set in the Z way of the
    matrix.

    Few reserved char:
    # (sharp): comment
    ! : to show the symbol and name of the files
    | (pipe char): which means the upper block continues vertically to
    here.
    _ (underline char): means that the previous structure continues
    here.
    - (dash): empty space

    Jul 08 2022
    Saeed
    """


super_str = struct.Structure()
files = super_str.files
param_fname = super_str.param_fname
output_fname = super_str.output
bs = mupdate.UpdateType(files, param_fname)
s = mstack.StackData(super_str.block, super_str.axis, bs)
wrt = lmpwrt.WriteLmp(s, output_fname)
wrt.write_lmp()
parm = parwrt.WriteParam(s, bs)
