import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from pysnptools.standardizer import Unit
from kernelreader import KernelReader
from kerneldata import KernelData

class SnpKernel(KernelReader):
    #!!!cmk update all comments
    """  This is a class hold SNP values in-memory along with related iid and sid information.
    It is created by calling the :meth:`.SnpReader.read` method on another :class:`.SnpReader`, for example, :class:`.Bed`.

    See :class:`.SnpReader` for details and examples.
    """
    def __init__(self, snpreader, standardizer=None, test=None, block_size=None): #!!!autodoc doesn't generate good doc for this constructor
        assert standardizer is not None, "'standardizer' must be provided"

        self.snpreader = snpreader
        self.standardizer = standardizer
        self.test = test or snpreader
        self.block_size = block_size

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
        return self.test.iid

    def __repr__(self):
        if self.snpreader is self.test:
            if self.block_size is None:
                s = "SnpKernel({0},standardizer={1})".format(self.snpreader,self.standardizer)
            else:
                s = "SnpKernel({0},standardizer={1},block_size={2})".format(self.snpreader,self.standardizer,self.block_size)
        else:
            if self.block_size is None:
                s = "SnpKernel({0},standardizer={1},test={2})".format(self.snpreader,self.standardizer,self.test)
            else:
                s = "SnpKernel({0},standardizer={1},test={2},block_size={3})".format(self.snpreader,self.standardizer,self.test,self.block_size)
        return s

    def copyinputs(self, copier):
        #Doesn't need run_once
        copier.input(self.snpreader)
        copier.input(self.standardizer)
        if self.snpreader is not self.test:
            copier.input(self.test)

    #!!!cmk0
    def _read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok):
        if row_index_or_none is col_index_or_none or np.array_equal(row_index_or_none,row_index_or_none):
            if row_index_or_none is None:
                return self.snpreader._read_kernel(self.standardizer,self.test,self.block_size,order, dtype, force_python_only, view_ok)
            else:
                #!!!CMK0
                raise NotImplementedError("Don't currently support reading non-square kernels from SnpKernels")
                return self.snpreader[row_index_or_none,:]._read_kernel(self.standardizer,self.test,self.block_size,order, dtype, force_python_only, view_ok)
        else:
            #!!!cmk0: Need to find the iids that are in either the cols or the rows. Standardize the snps with just those and then turn that square into the rectangle requested
            raise NotImplementedError("Don't currently support reading non-square kernels from SnpKernels")

    #!!!cmk0
    #!!!cmk be sure to document that any subsetting applies to the inter snpreader BEFORE standardization.
    def __getitem__(self, iid_indexer_and_snp_indexer):
        assert not isinstance(iid_indexer_and_snp_indexer, str), "Don't expect SnpKernel to be subsetted with a string" #LATER is this the best place for this test?
        try:
            iid0,iid1 = iid_indexer_and_snp_indexer
        except:
            iid0 = iid_indexer_and_snp_indexer
            iid1 = iid0

        #LATER is '==' and array_equal the right way to check for all possible advanced indexing, e.g. slices, etc.
        assert iid0 is iid1 or np.array_equal(iid0,iid1), "when selecting a subset of snps from a SnpKernel, the two snps lists must be the same" #LATER is this restriction good?

        return SnpKernel(self.snpreader[iid0,:],standardizer=self.standardizer)



if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
