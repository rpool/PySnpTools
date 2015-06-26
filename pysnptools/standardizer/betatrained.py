import numpy as np
import scipy as sp
import logging
import warnings
from pysnptools.standardizer import Standardizer

class BetaTrained(Standardizer):
    """description of class""" #!!!Cmk

    def __init__(self, a,b,stats,sid):
        self.a=a
        self.b=b
        self.stats=stats
        self.sid=sid

    def __repr__(self): 
        return "{0}(a={1},b={2},stats={3})".format(self.__class__.__name__,self.a,self.b,self.stats)

    #!!!cmk need doc
    def standardize(self, snps, block_size=None, force_python_only=False):
        if block_size is not None:
            warnings.warn("block_size is deprecated (and not needed, since standardization is in-place", DeprecationWarning)
        self._standardize_unit_and_beta(snps, is_beta=True, a=self.a, b=self.b, apply_in_place=True,use_stats=True,stats=self.stats,force_python_only=force_python_only)
        return snps
