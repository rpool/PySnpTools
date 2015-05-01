import numpy as np
import subprocess, sys
import os.path
from itertools import *
import pandas as pd
import logging
import time
import pysnptools.util as pstutil
import numbers

#!!why do the examples use ../tests/datasets instead of "examples"?
class PstReader(object):
    """The (abstract) base class for you to specify data and later read it.

    A PstReader is one of three things:

    #!!!cmk add in-memory examples
    * A :class:`.PstData` class that holds data in memory, typically after a read:

        >>> from pysnptools.pstreader import PstNpz
        >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz')
        >>> data1 = on_disk.read() #reads the values
        >>> type(data1.val) # The val property is an ndarray of values
        <type 'numpy.ndarray'>
        >>> print data1 # prints the specification of the in-memory information
        PstData(PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz'))
        >>> data1.row_count #prints the number of rows in this in-memory data
        300

    * A class such as :class:`.PstNpz` for you to specify data in file. For example,

        >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz')
        >>> print on_disk # prints specification for reading from file
        PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz')
        >>> on_disk.col_count # prints the number of columns (but doesn't read any column values)
        1015

    * A subset of any PstReader, specified with "[ *row_index* , *col_index* ]", to read just some values.

        >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz')
        >>> subset_on_disk = on_disk[[3,4],::2] # specification for a subset of the data on disk. No values are read yet.
        >>> print subset_on_disk.col_count # prints the number of columns in this subset (but still doesn't read any values)
        508
        >>> print subset_on_disk #prints a specification of 'subset_on_disk'
        PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz')[[3,4],::2]
        >>> data_subset = subset_on_disk.read() # efficiently reads the specified subset of values from the disk
        >>> print data_subset # prints the specification of the in-memory information
        PstData(PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz')[[3,4],::2])
        >>> int(data_subset.val.shape[0]),int(data_subset.val.shape[1]) # The dimensions of the ndarray of values
        (2, 508)

  
    Methods & Properties:

        Every PstReader, such as :class:`.PstNpz` and :class:`.PstData`, has these properties: :attr:`row`, :attr:`row_count`, :attr:`col`, :attr:`col_count`,
        :attr:`row_property`, :attr:`col_property` and these methods: :meth:`read`, :meth:`row_to_index`, :meth:`col_to_index`. See below for details.

        :class:`.PstData` is a PstReader so it supports the above properties and methods. In addition, it supports property :attr:`.PstData.val.
        See below for details.

    rows and cols:

        Example:

        >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz')
        >>> print on_disk.row[:3] # print the first three rows
        [['POP1' '0']
         ['POP1' '12']
         ['POP1' '44']]
        >>> print on_disk.col[:10] # print the first ten columns
        ['1_12' '1_34' '1_10' '1_35' '1_28' '1_25' '1_36' '1_39' '1_4' '1_13']
        >>> print on_disk.row_to_index([['POP1','44'],['POP1','12']]) #Find the indexes for two rows.
        [2 1]
        
    When PstData is Read:

        PstData can be enormous so we generally avoid reading it to the degree practical. Specifically,
        
        * Constructing and printing a PstReader causes no file reading. For example, these commands read no data:

            >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz') # Construct a PstNpz PstReader. No data is read.
            >>> print on_disk # Print the PstNpz PstReader specification. No data is read.
            PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz')
            >>> subset_on_disk = on_disk[[3,4],::2] # Construct a subsetting PstReader. No data is read.
            >>> print subset_on_disk # print the subset PstReader. No data is read.
            PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz')[[3,4],::2]

        * Properties and methods related to the rows and columns (to the degree practical) read only row and col data from the disk,
          not value data. Moreover, the row and col data is read from file only once. Consider these commands:

            >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz') # Construct a PstNpz PstReader. No data is read.
            >>> print on_disk.col[:10] # without reading any values data from disk, read the row and cik data from disk, cache it, and then print the first ten cols.
            ['1_12' '1_34' '1_10' '1_35' '1_28' '1_25' '1_36' '1_39' '1_4' '1_13']
            >>> print on_disk.col_to_index(['1_10','1_13']) #use the cached col information to find the indexes of '1_10' and '1_13'. (No data is read from disk.)
            [2 9]

        * The only methods that read values from file are :meth:`read (to the degree practical). For example:

            >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz') # Construct a PstNpz PstReader. No data is read.
            >>> data1 = on_disk.read() #read all the values from disk, creating a new PstData instance that keeps these values in memory
            >>> print data1.val[0,2] # print the value for the row with index 0 and the col with index 2. (No data is read from disk.)
            1.0

        * If you request the values for only a subset of the rows or columns, (to the degree practical) only that subset will be read from disk.
          for example:
          #!!!cmk check that these don't appear in this file: iid, col, PstNpz, SNP
            >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz')[[3,4],::2] # Construct a subsetting PstReader. No data is read.
            >>> data_subset = subset_on_disk.read() # from disk, read the values for the rows with index 3 and 4 AND cols with even numbered indexes.
            >>> print data_subset.val[0,2] # print the value with subset row index 0 and col index 2 (corresponding to row index 3 and col index 4 in the full data). No data is read from disk.
            2.0

    When PstData is Re-Read and Copied:

        Every time you call a PstReader's :meth:`read` method, the PstReader re-reads the value data and returns a new in-memory :class:`.PstData`
        (with :attr:`.PstData.val` property containing a new ndarray of the values).

        Here is an example of what not to do, because it causes all the SNP value data to be read twice.

            >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz') # Construct a PstNpz PstReader. No data is read.
            >>> # Not recommended because it inefficiently reads all the values twice.
            >>> print on_disk.read().val[0,2] # read all values into a new PstData, print a value
            1.0
            >>> print on_disk.read().val[0,3] # read all values (again) into a second new PstData, print a value
            2.0

        Here are two efficient alternatives. First, if all values can all fit in memory, read them once into a :class:`PstData` and then
        access that :class:`PstData` multiple times.

            >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz') # Construct a PstNpz PstReader. No data is read.
            >>> data1 = on_disk.read() # read all values into a new PstData
            >>> print data1.val[0,2] # print a value from data1's in-memory ndarray
            1.0
            >>> print data1.val[0,3] # print another value from data1's in-memory ndarray.
            2.0

        Second, if the value data is too large to fit in memory, use subsetting to read only the values of interest from disk.
       
            >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz') # Construct a PstNpz PstReader. No data is read.
            >>> print on_disk[0,2].read().val[0,0] #Define the subset of data and read only that subset from disk.
            1.0
            >>> print on_disk[0,3].read().val[0,0] #Define a second subset of data and read only that subset from disk.
            2.0

        Because the in-memory :class:`.PstData` class is a kind of PstReader, you may read from it, too.
        Doing so create a new :class:`.PstData` instance containing a copy of the SNP values in a new ndarray.

            >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz') # Construct a PstNpz PstReader. No data is read.
            >>> data1 = on_disk.read() # read all SNP values from disk into a new PstData
            >>> print data1.val is data1.val # Do the in-memory SNP values use the same memory as themselves? Yes
            True
            >>> data2 = data1.read() # copy all the SNP values into a new ndarray in a new PstData
            >>> print data2.val is data1.val # Do the two ndarrays of in-memory SNP values use the same memory?
            False


    Avoiding Unwanted ndarray Allocations

        You may want a subset of SNPs values from an in-memory :class:`PstData` and you may know that this subset and the original :class:`PstData`
        can safely share the memory of the ndarray of SNP values. For this case, the :meth:`read` has optional parameters called view_ok and order. If you override 
        the defaults of "view_ok=False,order='F'" with "view_ok=True,order='A', the :meth:`read` will, if practical, return a new 
        :class:`PstData` with a ndarray shares memory with the original ndarray.
        Use these parameters with care because any change to either ndarray (for example, via :meth:`.PstData.standardize`) will effect
        the others. Also keep in mind that :meth:`read` relies on ndarray's mechanisms to decide whether to actually
        share memory and so it may ignore your suggestion and allocate a new ndarray anyway.

            >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz') # Construct a PstNpz PstReader. No data is read.
            >>> data1 = on_disk.read() # read all data from disk into a PstData with a new ndarray
            >>> column01 = data1[:,0:1].read(view_ok=True,order='A') #create PstData with the data from just the first two SNPs. Sharing memory is OK. The memory may be laid out in any order (that is col-major and row-major are both OK).
            >>> import numpy as np
            >>> #print np.may_share_memory(data1.val, column01.val) # Do the two ndarray's share memory? They could (but currently they won't)
            >>> column201 = data1[:,[2,0,1]].read(view_ok=True,order='A') #create PstData with the data from three SNPs, permuted. Sharing memory is OK.
            >>> print np.may_share_memory(data1.val, column201.val) # Do the two ndarray's share memory? No, ndarray decided that this indexing was too complex for sharing.
            False

    Creating Subsetting PstReaders with Indexing

        You often don't want to read the SNP values for all rows and cols. You can use indexing to create a subsetting PstReader that
        will read only the SNP values of interest.

        PstReaders support the indexing formats supported by ndarray plus two generalizations. Here are examples of indexing with an array
        of indexes, with slicing, and with an array of Booleans.

            >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz') # Specify some data on disk in PstNpz format
            >>> subset_reader_1 = on_disk[[3,4],:] #index with an array of indexes
            >>> print subset_reader_1.row_count, subset_reader_1.col_count
            2 1015
            >>> data1 = subset_reader_1.read() # read just the two rows of interest from the disk
            >>> subset_reader_2 = on_disk[:,:0:-2] #index with a slice
            >>> print subset_reader_2.row_count, subset_reader_2.col_count
            300 507
            >>> boolindexes = [s.startswith('23_') for s in on_disk.col] # create a Boolean index of cols that start '23_'
            >>> subset_reader_3 = on_disk[:,boolindexes] #index with array of Booleans
            >>> print subset_reader_3.row_count, subset_reader_3.col_count
            300 24

        The first generalization over with ndarray offers is full indexing on both the row dimension and the col dimension, in other words,
        full multidimensional indexing. For example,

            >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz') # Specify some data on disk in PstNpz format
            >>> subset_reader_4= on_disk[[3,4],:0:-2] # index on two dimensions at once
            >>> print subset_reader_4.row_count, subset_reader_4.col_count
            2 507

        The second generalization is indexing on a single integer index.

            >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz') # Specify some data on disk in PstNpz format
            >>> subset_reader_5 = on_disk[5,:] #index with single integer
            >>> print subset_reader_5.row_count, subset_reader_5.col_count
            1 1015

        Indexing is also useful when you have SNP values in memory via a :class:`PstData` index and want to copy a subset of those values.
        While you could instead index directly on the `.PstData.val` ndarray, by indexing on the :class:`PstData` instance you
        also get row and cid information.

            >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz') # Specify some data on disk in PstNpz format
            >>> data1 = on_disk.read() # read all SNP values into memory
            >>> print data1.col[:10] # print the first 10 cols
            ['1_12' '1_34' '1_10' '1_35' '1_28' '1_25' '1_36' '1_39' '1_4' '1_13']
            >>> data_subset = data1[:,::2].read(view_ok=True,order='A') # create a copy or view with every other col
            >>> print data_subset.col[:10] # print the first 10 cols in the subset
            ['1_12' '1_10' '1_28' '1_36' '1_4' '1_11' '1_32' '1_9' '1_17' '1_18']


        You can apply indexing on top of indexing to specify subsets of subsets of data to read. In this example, 
        only the SNP values for every 16th col is actually read from the disk.

            >>> # These are just PstReaders, nothing is read from disk yet
            >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz') # Specify some data on disk in PstNpz format
            >>> half_reader = on_disk[:,::2] # a reader for half the cols
            >>> quarter_reader = half_reader[:,::2] # a reader for half of half the cols
            >>> sixteenth_reader = quarter_reader[:,::2][:,::2] # a reader for half of half of half of half the cols
            >>> print sixteenth_reader #Print the specification of this reader
            PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz')[:,::2][:,::2][:,::2][:,::2]
            >>> # Now we read from disk. Only values for one col in every 16 will be read.
            >>> data_sixteenth = sixteenth_reader.read()
            >>> print data_sixteenth.val[0,3]
            2.0

    The :meth:`read` Method
  
        By default the :meth:`read` returns a ndarray of scipy.float64 laid out in memory in F-contiguous order (row-index varies the fastest). You may, instead,
        ask for scipy.float32 or for C-contiguous order or any order. See :meth:`read` for details.
    """

    @property
    def row(self):
        """A ndarray of the rows. Each row is a ndarray of two strings (a family ID and a case ID) that identifies an individual.

        :rtype: ndarray (length :attr:`.row_count`) of ndarray (length 2) of strings

        This property (to the degree practical) reads only row and col data from the disk, not SNP value data. Moreover, the row and col data is read from file only once.

        :Example:

        >>> from pysnptools.pstreader import PstNpz
        >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz')
        >>> print on_disk.row[:3] # print the first three rows
        [['POP1' '0']
         ['POP1' '12']
         ['POP1' '44']]
        """
        raise NotImplementedError

    @property
    def row_count(self):
        """number of rows

        :rtype: integer

        This property (to the degree practical) reads only row and col data from the disk, not SNP value data. Moreover, the row and col data is read from file only once.
        """
        return len(self.row)

    @property
    def col(self):
        """A ndarray of the cols. Each col is a string that identifies a SNP.

        :rtype: ndarray (length :attr:`.col_count`) of strings

        This property (to the degree practical) reads only row and col data from the disk, not SNP value data. Moreover, the row and col data is read from file only once.

        :Example:

        >>> from pysnptools.pstreader import PstNpz
        >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz')
        >>> print on_disk.col[:10] # print the first ten cols
        ['1_12' '1_34' '1_10' '1_35' '1_28' '1_25' '1_36' '1_39' '1_4' '1_13']

        """
        raise NotImplementedError

    @property
    def col_count(self):
        """number of cols

        :rtype: integer

        This property (to the degree practical) reads only row and col data from the disk, not SNP value data. Moreover, the row and col data is read from file only once.

        """
        return len(self.col)

    #!!document that chr must not be X,Y,M only numbers (as per the PLINK PstNpz format)
    #!!Also what about telling the ref and alt allele? Also, what about tri and quad alleles, etc?
    @property
    def row_property(self):
        """A ndarray of the position information for each col. Each element is a ndarray of three scipy.numbers's (chromosome, genetic distance, basepair distance).

        :rtype: ndarray (length :attr:`.col_count`) of ndarray (length 3) of scipy.float64

        This property (to the degree practical) reads only row and col data from the disk, not SNP value data. Moreover, the row and col data is read from file only once.

        :Example:

        >>> from pysnptools.pstreader import PstNpz
        >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz')
        >>> # print the shape of the information for the first three rows. It's 3 x 0 because it there are three rows, but in this example there is 0 per row information.
        >>> print on_disk.row_property[:3].shape[0], on_disk.row_property[:3].shape[1]
        3 0
        """
        raise NotImplementedError

    #!!document that chr must not be X,Y,M only numbers (as per the PLINK PstNpz format)
    #!!Also what about telling the ref and alt allele? Also, what about tri and quad alleles, etc?
    @property
    def col_property(self):
        """A ndarray of the position information for each col. Each element is a ndarray of three scipy.numbers's (chromosome, genetic distance, basepair distance).

        :rtype: ndarray (length :attr:`.col_count`) of ndarray (length 3) of scipy.float64

        This property (to the degree practical) reads only row and col data from the disk, not SNP value data. Moreover, the row and col data is read from file only once.

        :Example:

        >>> from pysnptools.pstreader import PstNpz
        >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz')
        >>> print on_disk.col_property[:3] # print position information for the first three cols:
        [[ 1.          0.00800801  0.        ]
         [ 1.          0.023023    1.        ]
         [ 1.          0.0700701   4.        ]]
        """
        raise NotImplementedError


    def _read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok):
        raise NotImplementedError



    #!!check that views always return contiguous memory by default
    def read(self, order='F', dtype=np.float64, force_python_only=False, view_ok=False):
        """Reads the SNP values and returns a :class:`.PstData` (with :attr:`.PstData.val` property containing a new ndarray of the SNP values).


        :param order: {'F' (default), 'C', 'A'}, optional -- Specify the order of the ndarray. If order is 'F' (default),
            then the array will be in F-contiguous order (row-index varies the fastest).
            If order is 'C', then the returned array will be in C-contiguous order (col-index varies the fastest).
            If order is 'A', then the :attr:`.PstData.val`
            ndarray may be in any order (either C-, Fortran-contiguous, or even discontiguous).
        :type order: string or None

        :param dtype: {scipy.float64 (default), scipy.float32}, optional -- The data-type for the :attr:`.PstData.val` ndarray.
        :type dtype: data-type

        :param force_python_only: optional -- If False (default), may use outside library code. If True, requests that the read
            be done without outside library code.
        :type force_python_only: bool


        :param view_ok: optional -- If False (default), allocates new memory for the :attr:`.PstData.val`'s ndarray. If True,
            if practical and reading from a :class:`PstData`, will return a new 
            :class:`PstData` with a ndarray shares memory with the original :class:`PstData`.
            Typically, you'll also wish to use "order='A'" to increase the chance that sharing will be possible.
            Use these parameters with care because any change to either ndarray (for example, via :meth:`.PstData.standardize`) will effect
            the others. Also keep in mind that :meth:`read` relies on ndarray's mechanisms to decide whether to actually
            share memory and so it may ignore your suggestion and allocate a new ndarray anyway.
        :type view_ok: bool

        :rtype: :class:`.PstData`

        Calling the method again causes the SNP values to be re-read and creates a new in-memory :class:`.PstData` with a new ndarray of SNP values.

        If you request the values for only a subset of the rows or cols, (to the degree practical) only that subset will be read from disk.

        :Example:

        >>> from pysnptools.pstreader import PstNpz
        >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz') # Specify SNP data on disk
        >>> data1 = on_disk.read() # Read all the SNP data returning a PstData instance
        >>> print type(data1.val) # The PstData instance contains a ndarray of the data.
        <type 'numpy.ndarray'>
        >>> subset_data = on_disk[:,::2].read() # From the disk, read SNP values for every other col
        >>> print subset_data.val[0,0] # Print the first SNP value in the subset
        2.0
        >>> subsub_data = subset_data[:10,:].read(order='A',view_ok=True) # Create an in-memory subset of the subset with SNP values for the first ten rows. Share memory if practical.
        >>> import numpy as np
        >>> # print np.may_share_memory(subset_data.val, subsub_data.val) # Do the two ndarray's share memory? They could. Currently they won't.       
        """
        val = self._read(None, None, order, dtype, force_python_only, view_ok)
        from pstdata import PstData
        ret = PstData(self.row,self.col,self.row_property,self.col_property, val, str(self))
        return ret

    def row_to_index(self, list):
        """Takes a list of rows and returns a list of index numbers

        :param list: list of rows
        :type order: list of list of strings

        :rtype: ndarray of int
        
        This method (to the degree practical) reads only row and col data from the disk, not SNP value data. Moreover, the row and col data is read from file only once.

        :Example:

        >>> from pysnptools.pstreader import PstNpz
        >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz') # Specify SNP data on disk
        >>> print on_disk.row_to_index([['POP1','44'],['POP1','12']]) #Find the indexes for two rows.
        [2 1]
        """
        if not hasattr(self, "_row_to_index"):
            self._row_to_index = {}
            for index, item in enumerate(self.row):
                key = PstReader._makekey(item)
                if self._row_to_index.has_key(key) : raise Exception("Expect row to appear in data only once. ({0})".format(key))
                self._row_to_index[key] = index
        index = np.fromiter((self._row_to_index[PstReader._makekey(item1)] for item1 in list),np.int)
        return index

    def col_to_index(self, list):
        """Takes a list of cols and returns a list of index numbers

        :param list: list of cols
        :type list: list of strings

        :rtype: ndarray of int
        
        This method (to the degree practical) reads only row and col data from the disk, not SNP value data. Moreover, the row and col data is read from file only once.

        :Example:

        >>> from pysnptools.pstreader import PstNpz
        >>> on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz') # Specify SNP data on disk
        >>> print on_disk.col_to_index(['1_10','1_13']) #Find the indexes for two cols.
        [2 9]
        """
        if not hasattr(self, "_col_to_index"):
            self._col_to_index = {}
            for index, item in enumerate(self.col):
                key = PstReader._makekey(item)
                if self._col_to_index.has_key(key) : raise Exception("Expect snp to appear in data only once. ({0})".format(key))
                self._col_to_index[item] = index
        index = np.fromiter((self._col_to_index[PstReader._makekey(item1)] for item1 in list),np.int)
        return index

    @staticmethod
    def _makekey(item):
        try:
            hash(item)
            return item
        except:
            return tuple(item)

    def __getitem__(self, row_indexer_and_col_indexer):
        from _subset import _Subset
        row_indexer, col_indexer = row_indexer_and_col_indexer
        return _Subset(self, row_indexer, col_indexer)


    def copyinputs(self, copier):
        raise NotImplementedError

    @staticmethod
    def _is_all_slice(index_or_none):
        if index_or_none is None:
            return True
        return  isinstance(index_or_none,slice) and index_or_none == slice(None)

    @staticmethod
    def _make_sparray_or_slice(indexer):
        if indexer is None:
            return slice(None)

        if isinstance(indexer,np.ndarray):
            return PstReader._process_ndarray(indexer)

        if isinstance(indexer, slice):
            return indexer

        if np.isscalar(indexer):
            assert isinstance(indexer, numbers.Integral), "Expect scalar indexes to be integers"
            return np.array([indexer])

        return PstReader._process_ndarray(np.array(indexer))

    @staticmethod
    def _process_ndarray(indexer):
        if len(indexer)==0: # If it's very length the type is unreliable and unneeded.
            return np.zeros((0),dtype=np.integer)
        if indexer.dtype == bool:
            return np.arange(len(indexer),dtype=np.integer)[indexer]
        assert np.issubdtype(indexer.dtype, np.integer), "Indexer of unknown type"
        return indexer


    @staticmethod
    def _make_sparray_from_sparray_or_slice(count, indexer):
        if isinstance(indexer,slice):
            return apply(np.arange, indexer.indices(count))
        return indexer

    @staticmethod
    def _array_properties_are_ok(val, order, dtype):
        if val.dtype != dtype:
            return False
        if order is 'F':
            return val.flags['F_CONTIGUOUS']
        elif order is 'C':
            return val.flags['C_CONTIGUOUS']

        return True

    def _apply_sparray_or_slice_to_val(self, val, row_indexer_or_none, col_indexer_or_none, order, dtype, force_python_only):

        if (PstReader._is_all_slice(row_indexer_or_none) and PstReader._is_all_slice(col_indexer_or_none)  and not force_python_only and 
                (order == 'A' or (order == 'F' and val.flags['F_CONTIGUOUS']) or (order == 'C' and val.flags['C_CONTIGUOUS'])) and
                (dtype is None or  val.dtype == dtype)):
            return val, True

        row_indexer = PstReader._make_sparray_or_slice(row_indexer_or_none)
        col_indexer = PstReader._make_sparray_or_slice(col_indexer_or_none)
        if not force_python_only:
            row_index = PstReader._make_sparray_from_sparray_or_slice(self.row_count, row_indexer)
            col_index = PstReader._make_sparray_from_sparray_or_slice(self.col_count, col_indexer)
            sub_val = pstutil.sub_matrix(val, row_index, col_index, order=order, dtype=dtype)
            return sub_val, False


        if PstReader._is_all_slice(row_indexer) or PstReader._is_all_slice(col_indexer):
            sub_val = val[row_indexer, col_indexer] #!!is this faster than the C++?
        else: 
            row_index = PstReader._make_sparray_from_sparray_or_slice(self.row_count, row_indexer)
            col_index = PstReader._make_sparray_from_sparray_or_slice(self.col_count, col_indexer)
            #See http://stackoverflow.com/questions/21349133/numpy-array-integer-indexing-in-more-than-one-dimension
            sub_val = val[row_index.reshape(-1,1), col_index]

        assert len(sub_val.shape)==2, "Expect result of subsetting to be 2 dimensional"

        if not PstReader._array_properties_are_ok(sub_val, order, dtype):
            if order is None:
                order = "K"
            if dtype is None:
                dtype = sub_val.dtype
            sub_val = sub_val.astype(dtype, order, copy=True)

        shares_memory =  np.may_share_memory(val, sub_val)
        assert(PstReader._array_properties_are_ok(sub_val, order, dtype))
        return sub_val, shares_memory

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)


    #from pysnptools.pstreader import PstData
    #from pysnptools.snpreader import Bed
    #snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300').read()
    #data1 = PstData(row = snp_on_disk.row, #N x 2
    #                col = snp_on_disk.col, # S
    #                row_property = np.zeros([snp_on_disk.row_count,0]),  # N
    #                col_property = snp_on_disk.pos, # S x 3
    #                val = snp_on_disk.val # N x S
    #                )
    #print type(data1.val) # The val property is an ndarray of values
    ###<type 'numpy.ndarray'>
    ##print data1 # prints the specification of the in-memory information
    ###PstData()
    ##print data1.row_count #prints the number of rows in this in-memory data
    ###3

    #from pysnptools.pstreader import PstNpz
    #pstnpz_filename = '../tests/datasets/all_chr.maf0.001.N300.pst.npz'
    #PstNpz.write(data1,pstnpz_filename)
    #reader1 = PstNpz(pstnpz_filename)
    #data2 = reader1.read()

    #from pysnptools.pstreader import PstNpz
    #on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz')
    #print on_disk.row[:3] # print the first three rows
    ##[['POP1' '0']
    ## ['POP1' '12']
    ## ['POP1' '44']]
    #print on_disk.col[:10] # print the first ten columns
    ##['1_12' '1_34' '1_10' '1_35' '1_28' '1_25' '1_36' '1_39' '1_4' '1_13']
    #print on_disk.row_to_index([['POP1','44'],['POP1','12']]) #Find the indexes for two rows.
    ##[2 1]

    #from pysnptools.pstreader import PstNpz
    #on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz')
    #subset_on_disk = on_disk[[3,4],::2] # specification for a subset of the data on disk. No values are read yet.
    #print subset_on_disk.col_count # prints the number of columns in this subset (but still doesn't read any values)
    ##    508
    #print subset_on_disk #prints a specification of 'subset_on_disk'
    ##    PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz')[[3,4],::2]
    #data_subset = subset_on_disk.read() # efficiently reads the specified subset of values from the disk
    #print data_subset # prints the specification of the in-memory information
    ##    PstData(PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz')[[3,4],::2])
    #int(data_subset.val.shape[0]),int(data_subset.val.shape[1]) # The dimensions of the ndarray of values
    ##    (2, 508)

    #from pysnptools.pstreader import PstNpz
    #on_disk = PstNpz('../tests/datasets/all_chr.maf0.001.N300.pst.npz')
    #print on_disk.row_property[:3] # print position information for the first three cols:
    ##    [[ 1.          0.00800801  0.        ]
    ##     [ 1.          0.023023    1.        ]
    ##     [ 1.          0.0700701   4.        ]]


    #!!!cmk
    import doctest
    doctest.testmod()
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
    print "done"