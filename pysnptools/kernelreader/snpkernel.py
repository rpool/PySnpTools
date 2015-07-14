import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from pysnptools.standardizer import Unit
from kernelreader import KernelReader
from kerneldata import KernelData

class SnpKernel(KernelReader):
    '''
    A :class:`.KernelReader` that creates a kernel from a :class:`.SnpReader`. No SNP data will be read until
    the :meth:`SnpKernel.read` method is called. Use block_size to avoid ever reading all the SNP data into memory
    at once.

    See :class:`.KernelReader` for general examples of using KernelReaders.

    **Constructor:**
        :Parameters: * **snpreader** (:class:`SnpReader`) -- The SNP data
                     * **standardizer** (:class:`Standardizer`) -- How the SNP data should be standardized
                     * **test** (optional, :class:`SnpReader`) -- SNP data for test iids.
                     * **block_size** (optional, int) -- The number of SNPs to read at a time.

        If **test** is not given, then the kernel will be for the iids in **snpreader**.
        If **block_size** is not given, then all SNP data will be read at once.

        :Example:

        >>> from pysnptools.snpreader import Bed
        >>> from pysnptools.standardizer import Unit
        >>> snp_on_disk = Bed('../examples/toydata.bed')                    # A Bed file is specified, but nothing is read from disk
        >>> kernel_on_disk = SnpKernel(snp_on_disk, Unit(),block_size=500)  # A kernel is specified, but nothing is read from disk
        >>> print kernel_on_disk #Print the specification
        SnpKernel(Bed('../examples/toydata.bed'),standardizer=Unit(),block_size=500)
        >>> print kernel_on_disk.iid_count                                  # iid information is read from disk, but not SNP data
        500
        >>> kerneldata = kernel_on_disk.read().standardize()                # SNPs are read and Unit standardized, 500 at a time, to create a kernel, which is then standardized
        >>> print kerneldata.val[0,0]
        0.992306992842
    '''
    def __init__(self, snpreader, standardizer=None, test=None, block_size=None):
        assert standardizer is not None, "'standardizer' must be provided"

        self.snpreader = snpreader
        self.standardizer = standardizer
        self.test = test or snpreader
        self.block_size = block_size

    @property
    def row(self):
        return self.snpreader.iid

    @property
    def col(self):
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

    def _read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok):
        #Special case: If square and Identity, can push the subsetting into the SnpReader
        from pysnptools.standardizer import Identity as Stdizer_Identity
        if (isinstance(self.standardizer,Stdizer_Identity) and self.snpreader is self.test and 
            row_index_or_none is not None and col_index_or_none is not None and np.array_equal(row_index_or_none,col_index_or_none)):
            return self.snpreader[row_index_or_none,:]._read_kernel(self.standardizer,self.test,self.block_size,order, dtype, force_python_only, view_ok)
        else:
            #LATER: If it was often that case that we wanted to standardize on all the data, but then only return a slice of the result,
            #       that could be done with less memory by working in blocks but not tabulating for all the iids.
            whole = self.snpreader._read_kernel(self.standardizer,self.test,self.block_size,order, dtype, force_python_only, view_ok)
            val, shares_memory = self._apply_sparray_or_slice_to_val(whole, row_index_or_none, col_index_or_none, order, dtype, force_python_only)
            return val

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
