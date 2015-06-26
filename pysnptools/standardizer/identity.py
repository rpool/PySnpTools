import numpy as np
import scipy as sp
import logging
from pysnptools.standardizer import Standardizer

class Identity(Standardizer):
    """Do nothing to the data"""

    def __init__(self):
        pass

    def standardize(self, snps, block_size=None, force_python_only=False):
        if block_size is not None:
            warnings.warn("block_size is deprecated (and not needed, since standardization is in-place", DeprecationWarning)
        return snps

    #changes snpdata.val in place
    def _train_standardizer(self,snpdata,apply_in_place,force_python_only=False):
        return self


    def __repr__(self): 
        return "{0}()".format(self.__class__.__name__)


