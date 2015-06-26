import numpy as np
import scipy as sp
import logging
from pysnptools.standardizer import Standardizer

class BySidCount(Standardizer):
    """Standardize data by dividing every value by the number of SNPs"""
    def __init__(self, sid_count=None):
        self._sid_count = sid_count

    def standardize(self, snps, block_size=None, force_python_only=False): #!!!cmk0
        if block_size is not None:
            warnings.warn("block_size is deprecated (and not needed, since standardization is in-place", DeprecationWarning)

        sid_count = snps.shape[1] if self._sid_count is None else self._sid_count
        snps /= sid_count
        return snps

    def __repr__(self): 
            return "{0}()".format(self.__class__.__name__)
