import numpy as np
import scipy as sp
import logging
import warnings
from pysnptools.standardizer import Standardizer

class _CannotBeTrained(Standardizer):

    def __init__(self, name):
        self.name=name

    def __repr__(self): 
        return "{0}({1})".format(self.__class__.__name__,self.name)

    def standardize(self, snps, block_size=None, force_python_only=False):
        if block_size is not None:
            warnings.warn("block_size is deprecated (and not needed, since standardization is in-place", DeprecationWarning)
        raise Exception("Standardizer '{0}' cannot be trained")
