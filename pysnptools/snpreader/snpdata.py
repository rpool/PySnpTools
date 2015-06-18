import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from snpreader import SnpReader
from pysnptools.standardizer import Unit
from pysnptools.standardizer import Identity
from pysnptools.pstreader import PstData
import warnings

class SnpData(PstData,SnpReader):
    """  This is a class hold SNP values in-memory along with related iid and sid information.
    It is created by calling the :meth:`.SnpReader.read` method on another :class:`.SnpReader`, for example, :class:`.Bed`.

    See :class:`.SnpReader` for details and examples.
    """
    def __init__(self, iid, sid, val, pos=None, parent_string="", copyinputs_function=None): #!!!autodoc doesn't generate good doc for this constructor #!!!cmk5 should inits call the super inits to be sure everything is set?
        self._row = PstData._fixup_input(iid,empty_creator=lambda ignore:np.empty([0,2],dtype=str))
        self._col = PstData._fixup_input(sid,empty_creator=lambda ignore:np.empty([0],dtype=str))
        self._row_property = PstData._fixup_input(None,count=len(self._row),empty_creator=lambda count:np.empty([count,0],dtype=str))
        self._col_property = PstData._fixup_input(pos,count=len(self._col),empty_creator=lambda count:np.array([[np.nan, np.nan, np.nan]]*len(count)))
        self.val = PstData._fixup_input_val(val,row_count=len(self._row),col_count=len(self._col),empty_creator=lambda row_count,col_count:np.empty([row_count,col_count],dtype=np.float64))

        self._assert_iid_sid_pos()
        self._parent_string = parent_string
        self._std_string_list = []

    #LATER should there be a single warning if Unit() finds and imputes NaNs?
    def standardize(self, standardizer=Unit(), block_size=None, force_python_only=False):
        """Does in-place standardization of the in-memory
        SNP data. By default, it applies 'Unit' standardization, that is: the values for each SNP will have mean zero and standard deviation 1.0.
        NaN values are then filled with zero, the mean (consequently, if there are NaN values, the final standard deviation will not be zero.
        Note that, for efficiently, this method works in-place, actually changing values in the ndarray. Although it works in place, for convenience
        it also returns itself.

        :param standardizer: optional -- Specify standardization to be applied before the matrix multiply. 
             Any class from :mod:`pysnptools.standardizer` may be used. Some choices include :class:`.Unit` (default, makes values for each SNP have mean zero and
             standard deviation 1.0)and :class:`.Beta`.
        :type order: class from :mod:`pysnptools.standardizer`

        :param block_size: optional -- Default of None. None means to load all. Suggested number of sids to read into memory at a time.
        :type block_size: int or None

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
        self.val = standardizer.standardize(self.val, block_size=block_size, force_python_only=force_python_only)
        self._std_string_list.append(str(standardizer))
        return self

    def _read_kernel(self, standardizer, block_size=None, order='F', dtype=np.float64, force_python_only=False, view_ok=False):#!!!cmk0 respect all these inputs
        """
        See :meth:`.SnpReader.kernel` for details and examples.
        """
        if type(standardizer) is Identity:
            K = self.val.dot(self.val.T)
            return K
        else:
            K = SnpReader._read_kernel(self, standardizer, block_size=block_size, order=order, dtype=dtype, force_python_only=force_python_only,view_ok=view_ok)
            return K

    def __repr__(self):
        if self._parent_string == "":
            if len(self._std_string_list) > 0:
                s = "{0}({1})".format(self.__class__.__name__,",".join(self._std_string_list))
            else:
                s = "{0}()".format(self.__class__.__name__)
        else:
            if len(self._std_string_list) > 0:
                s = "{0}({1},{2})".format(self.__class__.__name__,self._parent_string,",".join(self._std_string_list))
            else:
                s = "{0}({1})".format(self.__class__.__name__,self._parent_string)
        return s


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
