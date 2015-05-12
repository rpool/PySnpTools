import logging
import scipy as sp
from kernelreader import KernelReader
from pysnptools.pstreader import PstHdf5
import warnings

#!! document the format

class KernelHdf5(PstHdf5,KernelReader):
    pass

