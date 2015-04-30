import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from kernelreader import KernelReader
from pstdata import PstData

class SnpKernel(KernelReader):
    #!!!cmk update all comments
    """  This is a class hold SNP values in-memory along with related iid and sid information.
    It is created by calling the :meth:`.SnpReader.read` method on another :class:`.SnpReader`, for example, :class:`.Bed`.

    See :class:`.SnpReader` for details and examples.
    """
    def __init__(self, snpreader, standardizer=Unit()): #!!!autodoc doesn't generate good doc for this constructor
        self.snpreader = snpreader
        self.standardizer = standardizer

    #!!!cmk does __repr__ do the right thing?
    #!!!cmk does copyinputs do the right thing? Probably not. Need to look at self.snpreader

    @property
    def assume_symmetric(self):
        return True

    @property
    def iid(self):
        """A ndarray of the iids.

        See :attr:`.SnpReader.iid` for details and examples.
        """
        return self.snpreader.iid

    @property
    def iid0(self):
        """A ndarray of the iid0s.

        See :attr:`.SnpReader.iid` for details and examples.
        """
        return self.snpreader.iid

    @property
    def iid1(self):
        """A ndarray of the iid0s.

        See :attr:`.SnpReader.iid` for details and examples.
        """
        return self.snpreader.iid

    #!!!cmk is this needed?
    # Most _read's support only indexlists or None, but this one supports Slices, too.
    _read_accepts_slices = None

    def _read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok):
        assert row_index_or_none == col_index_or_none, "!!!cmk fix up test and message"
        snpreader_subset = self.snpreader[row_index_or_none, col_index_or_none]
        val = snpreader_subset.kernel(self.standardizer) #!!!cmk what about order, dtype, and batch rows??
        return KernelData(iid=snpreader_subset.iid, val=val)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
