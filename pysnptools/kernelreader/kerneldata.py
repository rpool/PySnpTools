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
    def __init__(self, iid=None, iid0=None, iid1=None, val=None, parent_string="",copyinputs_function=None): #!!!autodoc doesn't generate good doc for this constructor
        assert val is not None, "'val' must not be None"
        assert iid is None ^ (iid0 is None and iid1 is None), "Either 'iid' or both 'iid0' 'iid1' must be provided."
        if iid is not None:
            iid0 = iid
            iid1 = iid

        self._row = iid0 if len(iid0)>0 else np.array([],dtype=str).reshape(0,2)
        self._col = iid1 if len(iid1)>0 else np.array([],dtype=str)

        self._assert_iid0_iid1()

        assert type(val) == np.ndarray, "expect SnpData's val to be a ndarray"
        self.val = val
        self._parent_string = parent_string
        self._std_string_list = []
        if copyinputs_function is not None:
            self.copyinputs = copyinputs_function

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
        assert self.assume_symmetric, "When 'iid' is used, iid0 must be the same as iid1"
        return self._iid0

    @property
    def iid0(self):
        """A ndarray of the iid0s.

        See :attr:`.SnpReader.iid` for details and examples.
        """
        return self._iid0

    @property
    def iid1(self):
        """A ndarray of the iid0s.

        See :attr:`.SnpReader.iid` for details and examples.
        """
        return self._iid1

    #!!!cmk is this needed?
    # Most _read's support only indexlists or None, but this one supports Slices, too.
    _read_accepts_slices = None



if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
