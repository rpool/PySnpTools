import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from kernelreader import KernelReader
from pysnptools.pstreader import PstData

class Identity(KernelReader):
    """  This is a class hold SNP values in-memory along with related iid and sid information.
    It is created by calling the :meth:`.SnpReader.read` method on another :class:`.SnpReader`, for example, :class:`.Bed`.

    See :class:`.SnpReader` for details and examples. #!!!cmkupdate
    """
    def __init__(self, iid): #!!!autodoc doesn't generate good doc for this constructor

        self._row = iid if len(iid)>0 else np.array([],dtype=str).reshape(0,2) #!!!cmk are these two lines right?

    """The in-memory SNP data. A numpy.ndarray with dimensions :attr:`.iid_count` x :attr:`.sid_count`.

    See :class:`.SnpReader` for details and examples.
    """

    def __repr__(self): 
        return "{0}()".format(self.__class__.__name__)


    @property
    def row(self):
        """A ndarray of the iid0s.

        See :attr:`.SnpReader.iid` for details and examples.
        """
        return self._row
    #!!!cmk don't we need col_propoerty and row_property, too?
    @property
    def col(self):
        """A ndarray of the iid0s.

        See :attr:`.SnpReader.iid` for details and examples.
        """
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
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
