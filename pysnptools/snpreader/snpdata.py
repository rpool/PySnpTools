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
    """A :class:`.SnpReader` for holding SNP values (or similar values) in-memory, along with related iid and sid information.
    It is usually created by calling the :meth:`.SnpReader.read` method on another :class:`.SnpReader`, for example, :class:`.Bed`. It can also be constructed.

    See :class:`.SnpReader` for details and examples.

    **Constructor:**
        :Parameters: * **iid** (an array of strings) -- The :attr:`.SnpReader.iid` information
                     * **sid** (an array of strings) -- The :attr:`.sid` information
                     * **val** (a 2-D array of floats) -- The SNP values
                     * **pos** (optional, an array of strings) -- The :attr:`.pos` information
                     * **parent_string** (optional, string) -- Information to be display about the origin of this data
                     * **copyinputs_function** (optional, function) -- *Used internally by optional clustering code*

        :Example:

        >>> from pysnptools.snpreader import SnpData
        >>> snpdata = SnpData(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'], val=[[0.,2.,0.],[0.,1.,2.]])
        >>> print snpdata.val[0,1], snpdata.iid_count, snpdata.sid_count
        2.0 2 3

    **Equality:**

        Two SnpData objects are equal if their four arrays (:attr:`.iid`, :attr:`.sid`, :attr:`.val`, and :attr:`.pos_property`) are 'array_equal'.
        (Their 'parent_string' does not need to be the same).

        :Example:

        >>> from pysnptools.snpreader import SnpData
        >>> snpdata1 = SnpData(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'], val=[[0.,2.,0.],[0.,1.,2.]], pos=[[0,0,0],[0,0,0],[0,0,0]])
        >>> snpdata2 = SnpData(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'], val=[[0.,2.,0.],[0.,1.,2.]], pos=[[0,0,0],[0,0,0],[0,0,0]])
        >>> print snpdata1 == snpdata2 #True, because all the arrays have the same values.
        True
        >>> print snpdata1.val is snpdata2.val #False, because the two arrays have different memory.
        False
        >>> snpdata3 = SnpData(iid=[['a','0'],['b','0']], sid=['snp334','snp349','snp921'], val=[[0.,2.,0.],[0.,1.,2.]], pos=[[0,0,0],[0,0,0],[0,0,0]])
        >>> snpdata4 = SnpData(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'], val=[[0.,2.,0.],[0.,1.,2.]], pos=[[0,0,0],[0,0,0],[0,0,0]])
        >>> print snpdata3 == snpdata4 #False, because the iids are different.
        False

    **Methods beyond** :class:`.SnpReader`
    """
    def __init__(self, iid, sid, val, pos=None, parent_string="", copyinputs_function=None):
        self._row = PstData._fixup_input(iid,empty_creator=lambda ignore:np.empty([0,2],dtype=str))
        self._col = PstData._fixup_input(sid,empty_creator=lambda ignore:np.empty([0],dtype=str))
        self._row_property = PstData._fixup_input(None,count=len(self._row),empty_creator=lambda count:np.empty([count,0],dtype=str))
        self._col_property = PstData._fixup_input(pos,count=len(self._col),empty_creator=lambda count:np.array([[np.nan, np.nan, np.nan]]*count))
        self.val = PstData._fixup_input_val(val,row_count=len(self._row),col_count=len(self._col),empty_creator=lambda row_count,col_count:np.empty([row_count,col_count],dtype=np.float64))

        self._assert_iid_sid_pos()
        self._parent_string = parent_string
        self._std_string_list = []

    val = None
    """The 2D NumPy array of floats that represents the values of the SNPs.

    >>> from pysnptools.snpreader import Bed
    >>> snpdata = Bed('../../tests/datasets/all_chr.maf0.001.N300')[:5,:].read() #read data for first 5 iids
    >>> print snpdata.val[4,100] #print one of the SNP values
    2.0
    """

    def train_standardizer(self, apply_in_place, standardizer=Unit(), force_python_only=False):
        """Trains a standardizer on the current in-memory SnpData that can be applied on other SNP data.

        :param apply_in_place: Tells if the current SnpData should be standardized.
        :type apply_in_place: bool

        :param standardizer: optional -- Specify standardization to be applied. 
             Any trainable :class:`.Standardizer` may be used. Some choices include :class:`.Unit` (default, makes values for each SNP have mean zero and
             standard deviation 1.0), :class:`.Beta`, and :class:`standardizer.Identity` (do nothing)
        :type standardizer: :class:`.Standardizer`

        :param force_python_only: optional -- If true, will use pure Python instead of faster C++ libraries.
        :type force_python_only: bool

        :rtype: :class:`.Standardizer` -- A trained standardizer

        >>> from pysnptools.snpreader import Bed
        >>> train = Bed('../../tests/datasets/all_chr.maf0.001.N300')[1:,:].read() # read SNP values for all but the first iid
        >>> train_stder = train.train_standardizer(apply_in_place=True) #Unit standardize and remember the mean and stddev of each SNP
        >>> test = Bed('../../tests/datasets/all_chr.maf0.001.N300')[0,:].read() # read SNP values for the first iid
        >>> test = test.standardize(train_stder) # Use the mean and stddev of the train data to unit standardize the test data.
        """
        if apply_in_place:
            self._std_string_list.append(str(standardizer))
        trained_standardizer = standardizer._train_standardizer(self, apply_in_place=apply_in_place, force_python_only=force_python_only)
        return trained_standardizer

    #LATER should there be a single warning if Unit() finds and imputes NaNs?
    def standardize(self, standardizer=Unit(), block_size=None, force_python_only=False):
        """Does in-place standardization of the in-memory
        SNP data. By default, it applies 'Unit' standardization, that is: the values for each SNP will have mean zero and standard deviation 1.0.
        NaN values are then filled with zero, the mean (consequently, if there are NaN values, the final standard deviation will not be zero.
        Note that, for efficiency, this method works in-place, actually changing values in the ndarray. Although it works in place, for convenience
        it also returns the SnpData.

        :param standardizer: optional -- Specify standardization to be applied. 
             Any :class:`.Standardizer` may be used. Some choices include :class:`.Unit` (default, makes values for each SNP have mean zero and
             standard deviation 1.0) and :class:`.Beta`.
        :type standardizer: :class:`.Standardizer`

        :param block_size: optional -- Deprecated.
        :type block_size: None

        :param force_python_only: optional -- If true, will use pure Python instead of faster C++ libraries.
        :type force_python_only: bool

        :rtype: :class:`.SnpData` (standardizes in place, but for convenience, returns 'self')

        >>> from pysnptools.snpreader import Bed
        >>> snp_on_disk = Bed('../../tests/datasets/all_chr.maf0.001.N300') # Specify some data on disk in Bed format
        >>> snpdata1 = snp_on_disk.read() # read all SNP values into memory
        >>> print snpdata1 # Prints the specification for this SnpData
        SnpData(Bed('../../tests/datasets/all_chr.maf0.001.N300'))
        >>> print snpdata1.val[0,0]
        2.0
        >>> snpdata1.standardize() # standardize changes the values in snpdata1.val and changes the specification.
        SnpData(Bed('../../tests/datasets/all_chr.maf0.001.N300'),Unit())
        >>> print snpdata1.val[0,0]
        0.229415733871
        >>> snpdata2 = snp_on_disk.read().standardize() # Read and standardize in one expression with only one ndarray allocated.
        >>> print snpdata2.val[0,0]
        0.229415733871
        """
        self._std_string_list.append(str(standardizer))
        _ = standardizer._train_standardizer(self, apply_in_place=True, force_python_only=force_python_only)
        return self

    def _read_kernel(train, standardizer, test=None, block_size=None, order='A', dtype=np.float64, force_python_only=False, view_ok=False):
        """
        See :meth:`.SnpReader.kernel` for details and examples.
        """
        from pysnptools.pstreader import PstReader

        if test is None:
            test = train

        #Just do a 'python' dot, if no standardization is needed and everything is the right type
        if isinstance(standardizer,Identity) and train.val.dtype == dtype and isinstance(test,SnpData) and test.val.dtype == dtype:
            if order == 'F': #numpy's 'dot' always returns 'C' order
                K = (test.val.dot(train.val.T)).T
            else:
                K = train.val.dot(test.val.T)
            assert PstReader._array_properties_are_ok(K,order,dtype), "internal error: K is not of the expected order or dtype"
            return K
        else: #Do things the more general SnpReader way.
            K = SnpReader._read_kernel(train, standardizer, test, block_size=block_size, order=order, dtype=dtype, force_python_only=force_python_only,view_ok=view_ok)
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
