import numpy as np
from pysnptools.standardizer import Standardizer

class DiagKtoN(Standardizer):
    """diag(K)=N standardization of the data"""
    def __init__(self, N):
        self._N = N

    def standardize(self, snps, block_size=None, force_python_only=False):#!!!cmk0
        if block_size is not None:
            warnings.warn("block_size is deprecated (and not needed, since standardization is in-place", DeprecationWarning)

        vec = snps.reshape(-1, order="A")
        
        # make sure no copy was made
        assert not vec.flags['OWNDATA']
        squared_sum = vec.dot(vec)
        factor = 1./np.sqrt(squared_sum / float(self._N))
        snps *= factor
                
        return snps

    def __repr__(self): 
        return "{0}()".format(self.__class__.__name__)
