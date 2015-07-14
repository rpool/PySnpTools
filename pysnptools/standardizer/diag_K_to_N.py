import numpy as np
from pysnptools.standardizer import Standardizer
import logging
import warnings

class DiagKtoN(Standardizer):
    '''
    A :class:`.Standardizer` that multiplies the SNP values by a factor such that a kernel
    constructed from the SNP values will have a diagonal that sums to iid_count. This class thus
    standardizes the kernel before it is even constructed.

    :Warning: This standardizer must run on all the SNPs values at once. If block_size is used, the wrong answer will result.

    >>> import numpy as np
    >>> from pysnptools.standardizer import DiagKtoN, Unit, Identity
    >>> from pysnptools.snpreader import Bed
    >>> snpdata1 = Bed('../../tests/datasets/all_chr.maf0.001.N300').read().standardize(Unit())
    >>> kernel1 = snpdata1.read_kernel(DiagKtoN(),block_size=None)
    >>> print np.diag(kernel1.val).sum()
    300.0
    '''
    """diag(K)=N standardization of the data"""
    def __init__(self,deprecated_iid_count=None):
        if deprecated_iid_count is not None:
            warnings.warn("'iid_count' is deprecated (and not needed, since can get iid_count from SNPs val's first dimension", DeprecationWarning)
        pass

    def standardize(self, snps, block_size=None, force_python_only=False):
        if block_size is not None:
            warnings.warn("block_size is deprecated (and not needed, since standardization is in-place", DeprecationWarning)

        vec = snps.reshape(-1, order="A")
        
        # make sure no copy was made
        assert not vec.flags['OWNDATA']
        squared_sum = vec.dot(vec)
        factor = 1./np.sqrt(squared_sum / float(snps.shape[0]))
        snps *= factor
                
        return snps

    def __repr__(self): 
        return "{0}()".format(self.__class__.__name__)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
        