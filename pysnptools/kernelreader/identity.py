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

    #_parent_string = ""
    #_std_string_list = []
    #def standardize(self):
    #    return self

    @property
    def iid(self):
        """A ndarray of the iids.

        See :attr:`.SnpReader.iid` for details and examples.
        """
        return self._row

    @property
    def iid0(self):
        """A ndarray of the iid0s.

        See :attr:`.SnpReader.iid` for details and examples.
        """
        return self._row

    @property
    def iid1(self):
        """A ndarray of the iid0s.

        See :attr:`.SnpReader.iid` for details and examples.
        """
        return self._row

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

    #!!!cmk test this
    def _read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok):
        #!!!cmk this code is not complete - if the row_index is not equal to the colum_index should slice the bit we want
        assert row_index_or_none == col_index_or_none, "!!!cmk fix up test and message"
        if row_index_or_none is None:
            return np.identity(self.row_count)
        else:
            return np.identity(len(row_index_or_none))

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
