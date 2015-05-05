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
        self._snpreader = snpreader
        self.standardizer = standardizer

    #!!!cmk does __repr__ do the right thing?
    #!!!cmk does copyinputs do the right thing? Probably not. Need to look at self.snpreader

    @property
    def row(self):
        """A ndarray of the iids.

        See :attr:`.SnpReader.iid` for details and examples.
        """
        return self._snpreader.iid

    @property
    def col(self):
        """A ndarray of the iid0s.

        See :attr:`.SnpReader.iid` for details and examples.
        """
        return self._snpreader.iid

    def __repr__(self):
        if isinstance(self.standardizer,Unit):
            s = "SnpKernel({0})".format(self._snpreader)
        else:
            s = "SnpKernel({0},standardizer={1})".format(self._snpreader,self.standardizer)
        return s

    def copyinputs(self, copier):
        #Doesn't need run_once
        copier.input(self._snpreader)
        copier.input(self.standardizer)


    #!!!cmk is this needed?
    # Most _read's support only indexlists or None, but this one supports Slices, too.
    _read_accepts_slices = None

    @property
    def sid_count(self):
        return self._snpreader.sid_count

    @property
    def sid(self):
        return self._snpreader.sid

    @property
    def pos(self):
        return self._snpreader.pos

    def read_snpdata(self, order='F', dtype=np.float64, force_python_only=False, view_ok=False):
        snpdata = self._snpreader.read(order=order, dtype=dtype, force_python_only=force_python_only, view_ok=view_ok)
        return snpdata.standardize(self.standardizer)

    def _read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok):
        #!!!cmk this code is not complete - if the row_index is not equal to the colum_index then should read the union of them, then compute the kernel then slice the bit we want
        assert np.array_equal(row_index_or_none, col_index_or_none), "!!!cmk fix up test and message"
        snpreader_subset = self._snpreader[row_index_or_none, :]
        #!!!cmk this doesn't seem right. The standardize (e.g. Unit) won't have the same effect after iids have been removed.
        val = snpreader_subset.kernel(self.standardizer) #!!!cmk what about order, dtype, and batch rows??
        return val

    #!!!cmk be sure to document that any subsetting applies to the inter snpreader BEFORE standardization.
    def __getitem__(self, iid_indexer_and_snp_indexer):
        assert not isinstance(iid_indexer_and_snp_indexer, str), "Don't expect SnpKernel to be subsetted with a string" #!!!is this the best place for this test?
        try:
            iid0,iid1 = iid_indexer_and_snp_indexer
        except:
            iid0 = iid_indexer_and_snp_indexer
            iid1 = iid0

        assert iid0 is iid1 or np.array_equal(iid0,iid1), "when selecting a subset of snps from a SnpKernel, the two snps lists must be the same" #!!!cmk is this restriction good?

        return SnpKernel(self._snpreader[iid0,:],standardizer=self.standardizer)



if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
