import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from kernelreader import KernelReader
from pysnptools.pstreader import PstData

class KernelData(KernelReader,PstData):
    """  This is a class hold SNP values in-memory along with related iid and sid information.
    It is created by calling the :meth:`.SnpReader.read` method on another :class:`.SnpReader`, for example, :class:`.Bed`.

    See :class:`.SnpReader` for details and examples.
    """
    def __init__(self, iid=None, iid0=None, iid1=None, val=None, parent_string=""): #!!!autodoc doesn't generate good doc for this constructor
        assert val is not None, "'val' must not be None"
        assert (iid is None) != (iid0 is None and iid1 is None), "Either 'iid' or both 'iid0' 'iid1' must be provided."
        if iid is not None:
            iid0 = iid
            iid1 = iid

        self._row = iid0 if len(iid0)>0 else np.array([],dtype=str).reshape(0,2) #!!!cmk are these two lines right?
        self._col = iid1 if len(iid1)>0 else np.array([],dtype=str)#!!!cmk are these two lines right? If empty, shouldn't they point to the SAME empty thing?

        assert type(val) == np.ndarray, "expect SnpData's val to be a ndarray"
        self.val = val
        self._parent_string = parent_string
        self._std_string_list = []

    val = None
    """The in-memory SNP data. A numpy.ndarray with dimensions :attr:`.iid_count` x :attr:`.sid_count`.

    See :class:`.SnpReader` for details and examples.
    """


    #!!!cmk does __repr__ do the right thing?
    #!!!cmk does copyinputs do the right thing?

    @property
    def iid(self):
        """A ndarray of the iids.

        See :attr:`.SnpReader.iid` for details and examples.
        """
        assert self.iid0 is self.iid1, "When 'iid' is used, iid0 must be the same as iid1"
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
        return self._col

    #!!!cmk is this needed?
    # Most _read's support only indexlists or None, but this one supports Slices, too.
    _read_accepts_slices = None

    def standardize(self):
        factor = self.iid_count / np.diag(self.val).sum()
        if abs(factor-1.0)>1e-15:
            self.val *= factor
        return self


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
