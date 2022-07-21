import os
import re
import pandas as pd
import numpy as np


class Doc:
    """prepare LAMMPS input file
    Preparing LAMMPS input file based on a base file (maybe!)
    The parameter input file will be in JSON format.
    The input file should contain:
        dt: timestep
        data: system data file
        f_field: forcefield file (parameter.lmp)
    """