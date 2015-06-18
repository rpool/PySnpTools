import logging
import scipy as sp
from snpreader import SnpReader
from pysnptools.pstreader import PstHdf5
import warnings

#!! document the format

class SnpHdf5(PstHdf5,SnpReader):
    pass

class Hdf5(SnpHdf5):
    warnings.warn("class 'Hdf5' is deprecated. Use the standard class 'SnpHdf5' instead", DeprecationWarning)
    pass

