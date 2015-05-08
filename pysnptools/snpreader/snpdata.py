import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from snpreader import SnpReader
from pysnptools.standardizer import Unit
from pysnptools.standardizer import Identity
from pysnptools.pstreader import PstData

class SnpData(PstData,SnpReader):
    """  This is a class hold SNP values in-memory along with related iid and sid information.
    It is created by calling the :meth:`.SnpReader.read` method on another :class:`.SnpReader`, for example, :class:`.Bed`.

    See :class:`.SnpReader` for details and examples.
    """
    def __init__(self, iid, sid, pos, val, parent_string="",copyinputs_function=None): #!!!autodoc doesn't generate good doc for this constructor #!!!cmk5 should inits call the super inits to be sure everything is set?
        self._row = iid if len(iid)>0 else np.array([],dtype=str).reshape(0,2) #!!!cmk4
        self._col = sid if len(sid)>0 else np.array([],dtype=str)
        self._col_property = pos if len(sid)>0 else np.array([],dtype=int).reshape(0,3) #!!!cmk4
        self._row_property = np.empty((len(iid),0))
        self._assert_iid_sid_pos()

        assert type(val) == np.ndarray, "expect SnpData's val to be a ndarray"
        self.val = val
        self._parent_string = parent_string
        self._std_string_list = []

    #!!! should there be a single warning if Unit() finds and imputes NaNs?
    def standardize(self, standardizer=Unit(), blocksize=None, force_python_only=False):
        """Does in-place standardization of the in-memory
        SNP data. By default, it applies 'Unit' standardization, that is: the values for each SNP will have mean zero and standard deviation 1.0.
        NaN values are then filled with zero, the mean (consequently, if there are NaN values, the final standard deviation will not be zero.
        Note that, for efficiently, this method works in-place, actually changing values in the ndarray. Although it works in place, for convenience
        it also returns itself.

        :param standardizer: optional -- Specify standardization to be applied before the matrix multiply. 
             Any class from :mod:`pysnptools.standardizer` may be used. Some choices include :class:`.Unit` (default, makes values for each SNP have mean zero and
             standard deviation 1.0)and :class:`.Beta`.
        :type order: class from :mod:`pysnptools.standardizer`

        :param blocksize: optional -- Default of None. None means to load all. Suggested number of sids to read into memory at a time.
        :type blocksize: int or None

        :rtype: :class:`.SnpData` (standardizes in place, but for convenience, returns 'self')

        >>> from pysnptools.snpreader import Bed
        >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300') # Specify some data on disk in Bed format
        >>> snpdata1 = snp_on_disk.read() # read all SNP values into memory
        >>> print snpdata1 # Prints the specification for this SnpData
        SnpData(Bed('../tests/datasets/all_chr.maf0.001.N300'))
        >>> print snpdata1.val[0,0]
        2.0
        >>> snpdata1.standardize() # standardize changes the values in snpdata1.val and changes the specification.
        SnpData(Bed('../tests/datasets/all_chr.maf0.001.N300'),Unit())
        >>> print snpdata1.val[0,0]
        0.229415733871
        >>> snpdata2 = snp_on_disk.read().standardize() # Read and standardize in one expression with only one ndarray allocated.
        >>> print snpdata2.val[0,0]
        0.229415733871
        """
        self.val = standardizer.standardize(self.val, blocksize=blocksize, force_python_only=force_python_only)
        self._std_string_list.append(str(standardizer))
        return self

    def kernel(self, standardizer, blocksize=10000, allowlowrank=False):
        """
            See :meth:`.SnpReader.kernel` for details and examples.
        """
        if type(standardizer) is Identity:
            K = self.val.dot(self.val.T)
            return K
        else:
            K = SnpReader.kernel(self, standardizer, blocksize=blocksize, allowlowrank=allowlowrank)
            return K

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
