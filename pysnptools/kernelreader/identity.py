import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from kernelreader import KernelReader
from pysnptools.pstreader import PstData

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
    '''
    def __init__(self, iid):

        if len(iid)>0:
            self._row = iid
        else:
            self._row = np.empty([0,2],dtype=str)

    def __repr__(self): 
        return "{0}()".format(self.__class__.__name__)


    @property
    def row(self):
        return self._row

    @property
    def col(self):
        return self._row

    def _read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok):
        if row_index_or_none is None and col_index_or_none is None:
            return np.identity(self.row_count,dtype=dtype)
        elif row_index_or_none is col_index_or_none or np.array_equal(row_index_or_none, col_index_or_none):
            return np.identity(len(row_index_or_none),dtype=dtype)
        else:
            #This is less efficient than it could be because it create a big identity matrix and then slices it.
            val, shares_memory = self._apply_sparray_or_slice_to_val(np.identity(self.row_count,dtype=dtype), row_index_or_none, col_index_or_none, order, dtype, force_python_only)
            return val


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
