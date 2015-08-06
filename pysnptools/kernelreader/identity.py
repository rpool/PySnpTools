import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from kernelreader import KernelReader
from pysnptools.pstreader import PstData
from pysnptools.pstreader import PstReader

class Identity(KernelReader):
    '''
    A :class:`.KernelReader` for that represents an identity matrix. No memory for the values is allocated until :meth:`Identity.read` is called.

    See :class:`.KernelReader` for general examples of using KernelReaders.

    **Constructor:**
        :Parameters: * **iid** (an array of strings) -- The :attr:`KernelReader.iid` information

        :Example:

        >>> from pysnptools.kernelreader import Identity
        >>> identity = Identity(iid=[['fam0','iid0'],['fam0','iid1']])
        >>> print identity.iid_count
        2
        >>> print identity.read().val
        [[ 1.  0.]
         [ 0.  1.]]

        >>> identity = Identity(iid=[['fam0','iid0'],['fam0','iid1'],['fam0','iid2']],test=[['fam0','iid1'],['fam0','iid3']])
        >>> print identity.iid0_count, identity.iid1_count
        3 2
        >>> print identity.read().val
        [[ 0.  0.]
         [ 1.  0.]
         [ 0.  0.]]

    '''
    def __init__(self, iid, test=None): #!!!cmk add docs and test for test

        if test is None:
            test = iid

        if len(iid)>0:
            self._row0 = iid
        else:
            self._row0 = _empty

        if len(test)>0:
            self._row1 = test
        else:
            self._row1 = _empty

    _empty = np.empty([0,2],dtype=str)

    def __repr__(self): 
        return "{0}()".format(self.__class__.__name__)


    @property
    def row(self):
        return self._row0

    @property
    def col(self):
        return self._row1

    def _read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok):
        if row_index_or_none is None and col_index_or_none is None and self._row0 is self._row1:
            return np.identity(self.row_count,dtype=dtype)
        elif (row_index_or_none is col_index_or_none or np.array_equal(row_index_or_none, col_index_or_none)) and self._row0 is self._row1:
            return np.identity(len(row_index_or_none),dtype=dtype)
        elif self._row0 is self._row1:
            #!!!cmk This is less efficient than it could be because it create a big identity matrix and then slices it.
            val, shares_memory = self._apply_sparray_or_slice_to_val(np.identity(self.row_count,dtype=dtype), row_index_or_none, col_index_or_none, order, dtype, force_python_only)
            return val
        else:
            #!!!cmk This is also less efficient than it could be because it create a big identity matrix and then slices it.

            #In about O(col_count + row_count) fill in zeros
            big = np.zeros([self.row_count,self.col_count],dtype=dtype)
            common = set([PstReader._makekey(x) for x in self.row]) & set([PstReader._makekey(x) for x in self.col])
            big[self.row_to_index(common),self.col_to_index(common)] = 1.0
            val, shares_memory = self._apply_sparray_or_slice_to_val(big, row_index_or_none, col_index_or_none, order, dtype, force_python_only)
            return val



if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
