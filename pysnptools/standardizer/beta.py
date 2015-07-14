import numpy as np
import scipy as sp
import logging
import warnings
from pysnptools.standardizer import Standardizer

class Beta(Standardizer):
    '''
    A :class:`.Standardizer` to beta standardize SNP data.

    See :class:`.Standardizer` for more information about standardization.

    **Constructor:**
        :Parameters: * **a** (*float*) -- The *a* parameter of the beta distribution
                     * **b** (*float*) -- The *b* parameter of the beta distribution

    >>> from pysnptools.standardizer import Beta
    >>> from pysnptools.snpreader import Bed
    >>> snpdata1 = Bed('../../tests/datasets/all_chr.maf0.001.N300').read().standardize(Beta(1,25))
    >>> print snpdata1.val[0,0]
    0.68080194805
    '''
    def __init__(self,a,b):
        self.a = a
        self.b = b

    def __repr__(self): 
        return "{0}(a={1},b={2})".format(self.__class__.__name__,self.a,self.b)

    #changes snpdata.val in place
    def _train_standardizer(self,snpdata,apply_in_place,force_python_only=False):
        from pysnptools.standardizer import BetaTrained
        stats=self._standardize_unit_and_beta(snpdata.val, is_beta=True, a=self.a, b=self.b, apply_in_place=apply_in_place,use_stats=False,stats=None,force_python_only=force_python_only)
        return BetaTrained(self.a,self.b,stats)

    def standardize(self, snps, block_size=None, force_python_only=False):
        if block_size is not None:
            warnings.warn("block_size is deprecated (and not needed, since standardization is in-place", DeprecationWarning)
        self._standardize_unit_and_beta(snps, is_beta=True, a=self.a, b=self.b, apply_in_place=True, use_stats=False,stats=None,force_python_only=force_python_only)
        return snps

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
