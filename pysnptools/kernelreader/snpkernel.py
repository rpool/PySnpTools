import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from pysnptools.standardizer import Unit
from kernelreader import KernelReader
from kerneldata import KernelData
#from pstdata import PstData !!!cmk

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
    def row(self):
        """A ndarray of the iids.

        See :attr:`.SnpReader.iid` for details and examples.
        """
        return self.snpreader.iid

    @property
    def col(self):
        """A ndarray of the iid0s.

        See :attr:`.SnpReader.iid` for details and examples.
        """
        return self.snpreader.iid

    def __repr__(self):
        if isinstance(self.standardizer,Unit):
            s = "SnpKernel({0})".format(self.snpreader)
        else:
            s = "SnpKernel({0},standardizer={1})".format(self.snpreader,self.standardizer)
        return s


    #!!!cmk is this needed?
    # Most _read's support only indexlists or None, but this one supports Slices, too.
    _read_accepts_slices = None

    def _read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok):
        #!!!cmk this code is not complete - if the row_index is not equal to the colum_index then should read the union of them, then compute the kernel then slice the bit we want
        try:
            assert np.array_equal(row_index_or_none, col_index_or_none), "!!!cmk fix up test and message"
        except:
            print "!!!cmk"
            raise Exception()
        snpreader_subset = self.snpreader[row_index_or_none, :]
        val = snpreader_subset.kernel(self.standardizer) #!!!cmk what about order, dtype, and batch rows??
        return val


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
