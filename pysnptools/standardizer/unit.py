import numpy as np
import scipy as sp
import logging
from pysnptools.standardizer import Standardizer
import warnings

class Unit(Standardizer):
    """A :class:`.Standardizer` to unit standardize SNP data. For each sid, the mean of the values is zero with standard deviation 1.
    NaN values are then filled with zero, the mean (consequently, if there are NaN values, the final standard deviation will not be zero.

    See :class:`.Standardizer` for more information about standardization.

    >>> from pysnptools.standardizer import Unit
    >>> from pysnptools.snpreader import Bed
    >>> snpdata1 = Bed('../../tests/datasets/all_chr.maf0.001.N300').read().standardize(Unit())
    >>> print snpdata1.val[0,0]
    0.229415733871
    """
    def __init__(self):
        pass

    def __repr__(self): 
        return "{0}()".format(self.__class__.__name__)

    #changes snpdata.val in place
    def _train_standardizer(self,snpdata,apply_in_place,force_python_only=False):
        from pysnptools.standardizer import UnitTrained
        stats=self._standardize_unit_and_beta(snpdata.val, is_beta=False, a=np.nan, b=np.nan, apply_in_place=apply_in_place,use_stats=False,stats=None, force_python_only=force_python_only)
        return UnitTrained(stats)

    def standardize(self, snps, block_size=None, force_python_only=False):
        if block_size is not None:
            warnings.warn("block_size is deprecated (and not needed, since standardization is in-place", DeprecationWarning)
        self._standardize_unit_and_beta(snps, is_beta=False, a=np.nan, b=np.nan, apply_in_place=True,use_stats=False,stats=None,force_python_only=force_python_only,)
        return snps

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
