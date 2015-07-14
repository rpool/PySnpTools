import numpy as np
import scipy as sp
import logging
import warnings
from pysnptools.standardizer import Standardizer

class UnitTrained(Standardizer):
    """A :class:`.Standardizer` to unit standardize one set of SNP data based on the mean and stddev of another set of SNP data
    NaN values are then filled with zero.

    See :class:`.Standardizer` for more information about standardization.

    **Constructor:**
        :Parameters: * **stats** (*ndarray of float*) -- The mean and stddev of each sid

    >>> from pysnptools.standardizer import Unit
    >>> from pysnptools.snpreader import Bed
    >>> train = Bed('../../tests/datasets/all_chr.maf0.001.N300')[1:,:].read() # read SNP values for all but the first iid
    >>> unittrained = train.train_standardizer(apply_in_place=True,standardizer=Unit()) #Unit standardize and remember the mean and stddev of each sid
    >>> print unittrained.stats[:5,:] #Print the means and stddev of the first five sids
    [[ 1.94983278  0.21828988]
     [ 1.96989967  0.17086341]
     [ 1.84280936  0.39057474]
     [ 1.99665552  0.0577347 ]
     [ 1.97658863  0.15120608]]
    >>> test = Bed('../../tests/datasets/all_chr.maf0.001.N300')[0,:].read() # read SNP values for the first iid
    >>> test = test.standardize(unittrained) # Use the mean and stddev of the train data to unit standardize the test data.
    >>> print test.val[0,0]
    0.229819279888
    """

    #!!might want to add code so that can check that sids are in the same order for both test and train
    def __init__(self, stats):
        self.stats=stats

    def __repr__(self): 
        return "{0}(stats={1})".format(self.__class__.__name__,self.stats)

    def standardize(self, snps, block_size=None, force_python_only=False):
        if block_size is not None:
            warnings.warn("block_size is deprecated (and not needed, since standardization is in-place", DeprecationWarning)
        _ = self._standardize_unit_and_beta(snps, is_beta=False, a=np.nan, b=np.nan, apply_in_place=True,use_stats=True,stats=self.stats,force_python_only=force_python_only)
        return snps

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
