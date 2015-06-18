import numpy as np
import subprocess, sys
import os.path
from itertools import *
import pandas as pd
import logging
import time
import pysnptools.util as pstutil
from pysnptools.pstreader import PstReader

#!!why do the examples use ../tests/datasets instead of "examples"?
class KernelReader(PstReader):
    """The (abstract) base class for you to specify SNP data and later read it.

    A SnpReader is one of three things:

    * A class such as :class:`.Bed` for you to specify data in file. For example,

        >>> from pysnptools.snpreader import Bed
        >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300')
        >>> print snp_on_disk # prints specification for reading from file
        Bed('../tests/datasets/all_chr.maf0.001.N300')
        >>> snp_on_disk.sid_count # prints the number of SNPS (but doesn't read any SNP values)
        1015

    * A :class:`.SnpData` class that holds SNP data in memory, typically after a read:

        >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300')
        >>> snpdata1 = snp_on_disk.read() #reads the SNP values
        >>> type(snpdata1.val) # The val property is an ndarray of SNP values
        <type 'numpy.ndarray'>
        >>> print snpdata1 # prints the specification of the in-memory SNP information
        SnpData(Bed('../tests/datasets/all_chr.maf0.001.N300'))
        >>> snpdata1.iid_count #prints the number of iids (number of individuals) in this in-memory data
        300


    * A subset of any SnpReader, specified with "[ *iid_index* , *sid_index* ]", to read only some SNP values.

        >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300')
        >>> subset_on_disk = snp_on_disk[[3,4],::2] # specification for a subset of the data on disk. No SNP values are read yet.
        >>> print subset_on_disk.sid_count # prints the number of sids in this subset (but still doesn't read any SNP values)
        508
        >>> print subset_on_disk #prints a specification of 'subset_on_disk'
        Bed('../tests/datasets/all_chr.maf0.001.N300')[[3,4],::2]
        >>> snpdata_subset = subset_on_disk.read() # efficiently reads the specified subset of values from the disk
        >>> print snpdata_subset # prints the specification of the in-memory SNP information
        SnpData(Bed('../tests/datasets/all_chr.maf0.001.N300')[[3,4],::2])
        >>> int(snpdata_subset.val.shape[0]),int(snpdata_subset.val.shape[1]) # The dimensions of the ndarray of SNP values
        (2, 508)

  
    Methods & Properties:

        Every SnpReader, such as :class:`.Bed` and :class:`.SnpData`, has these properties: :attr:`iid`, :attr:`iid_count`, :attr:`sid`, :attr:`sid_count`,
        :attr:`pos` and these methods: :meth:`read`, :meth:`iid_to_index`, :meth:`sid_to_index`, :meth:`kernel`. See below for details.

        :class:`.SnpData` is a SnpReader so it supports the above properties and method. In addition, it supports property :attr:`.SnpData.val` and method :meth:`.SnpData.standardize`.
        See below for details.

    iids and sids:

        Individual are identified with an iid, which is a ndarray of two strings: a family ID and a case ID. SNP locations 
        are identified with sid string. For example:

        >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300')
        >>> print snp_on_disk.iid[:3] # print the first three iids
        [['POP1' '0']
         ['POP1' '12']
         ['POP1' '44']]
        >>> print snp_on_disk.sid[:10] # print the first ten sids
        ['1_12' '1_34' '1_10' '1_35' '1_28' '1_25' '1_36' '1_39' '1_4' '1_13']
        >>> print snp_on_disk.iid_to_index([['POP1','44'],['POP1','12']]) #Find the indexes for two iids.
        [2 1]
        
    When Data is Read:

        SNP data can be enormous so generally avoid reading it to the degree practical. Specifically,
        
        * Constructing and printing a SnpReader causes no file reading. For example, these commands read no data:

            >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300') # Construct a Bed SnpReader. No data is read.
            >>> print snp_on_disk # Print the Bed SnpReader specification. No data is read.
            Bed('../tests/datasets/all_chr.maf0.001.N300')
            >>> subset_on_disk = snp_on_disk[[3,4],::2] # Construct a subsetting SnpReader. No data is read.
            >>> print subset_on_disk # print the subset SnpReader. No data is read.
            Bed('../tests/datasets/all_chr.maf0.001.N300')[[3,4],::2]

        * Properties and methods related to the iids and sids (to the degree practical) read only iid and sid data from the disk,
          not SNP value data. Moreover, the iid and sid data is read from file only once. Consider these commands:

            >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300') # Construct a Bed SnpReader. No data is read.
            >>> print snp_on_disk.sid[:10] # without reading any SNP values data from disk, read the sid and iid data from disk, cache it, and then print the first ten sids.
            ['1_12' '1_34' '1_10' '1_35' '1_28' '1_25' '1_36' '1_39' '1_4' '1_13']
            >>> print snp_on_disk.sid_to_index(['1_10','1_13']) #use the cached sid information to find the indexes of '1_10' and '1_13'. (No data is read from disk.)
            [2 9]

        * The only methods that read SNP values from file are :meth:`read` and :meth:`kernel` (to the degree practical). For example:

            >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300') # Construct a Bed SnpReader. No data is read.
            >>> snpdata1 = snp_on_disk.read() #read all the SNP values from disk, creating a new SnpData instance that keeps these values in memory
            >>> print snpdata1.val[0,2] # print the SNP value for the iid with index 0 and the sid with index 2. (No data is read from disk.)
            1.0

        * If you request the values for only a subset of the sids or iids, (to the degree practical) only that subset will be read from disk.
          for example:

            >>> subset_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300')[[3,4],::2] # Construct a subsetting SnpReader. No data is read.
            >>> snpdata_subset = subset_on_disk.read() # from disk, read the SNP values for the iids with index 3 and 4 AND sids with even numbered indexes.
            >>> print snpdata_subset.val[0,2] # print the SNP value with subset iid index 0 and sid index 2 (corresponding to iid index 3 and sid index 4 in the full data). No data is read from disk.
            2.0

    When Data is Re-Read and Copied:

        Every time you call a SnpReader's :meth:`read` method, the SnpReader re-reads the SNP value data and returns a new in-memory :class:`.SnpData`
        (with :attr:`.SnpData.val` property containing a new ndarray of the SNP values). Likewise, when you call the :meth:`kernel` method, the SnpReader re-reads
        the data and returns a new kernel ndarray.

        Here is an example of what not to do, because it causes all the SNP value data to be read twice.

            >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300') # Construct a Bed SnpReader. No data is read.
            >>> # Not recommended because it inefficiently reads all the SNP values twice.
            >>> print snp_on_disk.read().val[0,2] # read all values into a new SnpData, print a SNP value
            1.0
            >>> print snp_on_disk.read().val[0,3] # read all values (again) into a second new SnpData, print a SNP value
            2.0

        Here are two efficient alternatives. First, if all SNP values can all fit in memory, read them once into a :class:`SnpData` and then
        access that :class:`SnpData` multiple times.

            >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300') # Construct a Bed SnpReader. No data is read.
            >>> snpdata1 = snp_on_disk.read() # read all values into a new SnpData
            >>> print snpdata1.val[0,2] # print a SNP value from snpdata1's in-memory ndarray
            1.0
            >>> print snpdata1.val[0,3] # print another SNP value from snpdata1's in-memory ndarray.
            2.0

        Second, if the SNP value data is too large to fit in memory, use subsetting to read only the SNP values of interest from disk.
       
            >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300') # Construct a Bed SnpReader. No data is read.
            >>> print snp_on_disk[0,2].read().val[0,0] #Define the subset of data and read only that subset from disk.
            1.0
            >>> print snp_on_disk[0,3].read().val[0,0] #Define a second subset of data and read only that subset from disk.
            2.0

        Because the in-memory :class:`.SnpData` class is a kind of SnpReader, you may read from it, too.
        Doing so create a new :class:`.SnpData` instance containing a copy of the SNP values in a new ndarray.

            >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300') # Construct a Bed SnpReader. No data is read.
            >>> snpdata1 = snp_on_disk.read() # read all SNP values from disk into a new SnpData
            >>> print snpdata1.val is snpdata1.val # Do the in-memory SNP values use the same memory as themselves? Yes
            True
            >>> snpdata2 = snpdata1.read() # copy all the SNP values into a new ndarray in a new SnpData
            >>> print snpdata2.val is snpdata1.val # Do the two ndarrays of in-memory SNP values use the same memory?
            False


    Avoiding Unwanted ndarray Allocations

        You may want a subset of SNPs values from an in-memory :class:`SnpData` and you may know that this subset and the original :class:`SnpData`
        can safely share the memory of the ndarray of SNP values. For this case, the :meth:`read` has optional parameters called view_ok and order. If you override 
        the defaults of "view_ok=False,order='F'" with "view_ok=True,order='A', the :meth:`read` will, if practical, return a new 
        :class:`SnpData` with a ndarray shares memory with the original ndarray.
        Use these parameters with care because any change to either ndarray (for example, via :meth:`.SnpData.standardize`) will effect
        the others. Also keep in mind that :meth:`read` relies on ndarray's mechanisms to decide whether to actually
        share memory and so it may ignore your suggestion and allocate a new ndarray anyway.

            >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300') # Construct a Bed SnpReader. No data is read.
            >>> snpdata1 = snp_on_disk.read() # read all data from disk into a SnpData with a new ndarray
            >>> column01 = snpdata1[:,0:1].read(view_ok=True,order='A') #create SnpData with the data from just the first two SNPs. Sharing memory is OK. The memory may be laid out in any order (that is sid-major and iid-major are both OK).
            >>> import numpy as np
            >>> #print np.may_share_memory(snpdata1.val, column01.val) # Do the two ndarray's share memory? They could (but currently they won't)
            >>> column201 = snpdata1[:,[2,0,1]].read(view_ok=True,order='A') #create SnpData with the data from three SNPs, permuted. Sharing memory is OK.
            >>> print np.may_share_memory(snpdata1.val, column201.val) # Do the two ndarray's share memory? No, ndarray decided that this indexing was too complex for sharing.
            False

    Creating Subsetting SnpReaders with Indexing

        You often don't want to read the SNP values for all iids and sids. You can use indexing to create a subsetting SnpReader that
        will read only the SNP values of interest.

        SnpReaders support the indexing formats supported by ndarray plus two generalizations. Here are examples of indexing with an array
        of indexes, with slicing, and with an array of Booleans.

            >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300') # Specify some data on disk in Bed format
            >>> subset_snpreader_1 = snp_on_disk[[3,4],:] #index with an array of indexes
            >>> print subset_snpreader_1.iid_count, subset_snpreader_1.sid_count
            2 1015
            >>> snpdata1 = subset_snpreader_1.read() # read just the two rows of interest from the disk
            >>> subset_snpreader_2 = snp_on_disk[:,:0:-2] #index with a slice
            >>> print subset_snpreader_2.iid_count, subset_snpreader_2.sid_count
            300 507
            >>> boolindexes = [s.startswith('23_') for s in snp_on_disk.sid] # create a Boolean index of sids that start '23_'
            >>> subset_snpreader_3 = snp_on_disk[:,boolindexes] #index with array of Booleans
            >>> print subset_snpreader_3.iid_count, subset_snpreader_3.sid_count
            300 24

        The first generalization over with ndarray offers is full indexing on both the iid dimension and the sid dimension, in other words,
        full multidimensional indexing. For example,

            >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300') # Specify some data on disk in Bed format
            >>> subset_snpreader_4= snp_on_disk[[3,4],:0:-2] # index on two dimensions at once
            >>> print subset_snpreader_4.iid_count, subset_snpreader_4.sid_count
            2 507

        The second generalization is indexing on a single integer index.

            >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300') # Specify some data on disk in Bed format
            >>> subset_snpreader_5 = snp_on_disk[5,:] #index with single integer
            >>> print subset_snpreader_5.iid_count, subset_snpreader_5.sid_count
            1 1015

        Indexing is also useful when you have SNP values in memory via a :class:`SnpData` index and want to copy a subset of those values.
        While you could instead index directly on the `.SnpData.val` ndarray, by indexing on the :class:`SnpData` instance you
        also get iid and cid information.

            >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300') # Specify some data on disk in Bed format
            >>> snpdata1 = snp_on_disk.read() # read all SNP values into memory
            >>> print snpdata1.sid[:10] # print the first 10 sids
            ['1_12' '1_34' '1_10' '1_35' '1_28' '1_25' '1_36' '1_39' '1_4' '1_13']
            >>> snpdata_subset = snpdata1[:,::2].read(view_ok=True,order='A') # create a copy or view with every other sid
            >>> print snpdata_subset.sid[:10] # print the first 10 sids in the subset
            ['1_12' '1_10' '1_28' '1_36' '1_4' '1_11' '1_32' '1_9' '1_17' '1_18']


        You can apply indexing on top of indexing to specify subsets of subsets of data to read. In this example, 
        only the SNP values for every 16th sid is actually read from the disk.

            >>> # These are just SnpReaders, nothing is read from disk yet
            >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300') # Specify some data on disk in Bed format
            >>> half_snpreader = snp_on_disk[:,::2] # a reader for half the sids
            >>> quarter_snpreader = half_snpreader[:,::2] # a reader for half of half the sids
            >>> sixteenth_snpreader = quarter_snpreader[:,::2][:,::2] # a reader for half of half of half of half the sids
            >>> print sixteenth_snpreader #Print the specification of this reader
            Bed('../tests/datasets/all_chr.maf0.001.N300')[:,::2][:,::2][:,::2][:,::2]
            >>> # Now we read from disk. Only values for one sid in every 16 will be read.
            >>> snpdata_sixteenth = sixteenth_snpreader.read()
            >>> print snpdata_sixteenth.val[0,3]
            2.0

    The :meth:`read` Method
  
        By default the :meth:`read` returns a ndarray of scipy.float64 laid out in memory in F-contiguous order (iid-index varies the fastest). You may, instead,
        ask for scipy.float32 or for C-contiguous order or any order. See :meth:`read` for details.

    The :meth:`.SnpData.standardize` Method
        The :meth:`.SnpData.standardize` method, available only on :class:`.SnpData`, does in-place standardization of the in-memory
        SNP data. By default, it applies 'Unit' standardization, that is: the values for each SNP will have mean zero and standard deviation 1.0.
        NaN values are then filled with zero, the mean (consequently, if there are NaN values, the final standard deviation will not be zero.
        Note that, for efficiently, this method works in-place, actually changing values in the ndarray. Although it works in place, for convenience
        it also returns itself. See :meth:`.SnpData.standardize` for options and details.

            >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300') # Specify some data on disk in Bed format
            >>> snpdata1 = snp_on_disk.read() # read all SNP values into memory
            >>> print snpdata1 # Prints the specification for this SnpData
            SnpData(Bed('../tests/datasets/all_chr.maf0.001.N300'))
            >>> print snpdata1.val[0,0]
            2.0
            >>> snpdata1.standardize() # standardize changes the values in snpdata1.val and changes the specification.
            SnpData(Bed('../tests/datasets/all_chr.maf0.001.N300'),Unit())
            >>> print snpdata1.val[0,0]
            0.229415733871
            >>> snpdata2 = snp_on_disk.read().standardize() # Read and standardize in one expression with only one ndarray allocated.
            >>> print snpdata2.val[0,0]
            0.229415733871
   
    The :meth:`kernel` Method

        The :meth:`kernel` method, available on any SnpReader, returns a ndarray of size iid_count x iid_count. The returned array has the value
        of the (optionally standardized) SNP values transposed and then multiplied with themselves. When applied to an read-from-disk SnpReader, such as :class:`.Bed`,
        the method can save memory by reading (and standardizing) the data in blocks. See :meth:`kernel` for details.


    """

    @property
    def iid(self):
        """A ndarray of the iids. Each iid is a ndarray of two strings (a family ID and a case ID) that identifies an individual.

        :rtype: ndarray (length :attr:`.iid_count`) of ndarray (length 2) of strings

        This property (to the degree practical) reads only iid and sid data from the disk, not SNP value data. Moreover, the iid and sid data is read from file only once.

        :Example:

        >>> from pysnptools.snpreader import Bed
        >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300')
        >>> print snp_on_disk.iid[:3] # print the first three iids
        [['POP1' '0']
         ['POP1' '12']
         ['POP1' '44']]
        """
        assert self.iid0 is self.iid1, "When 'iid' is used, iid0 must be the same as iid1"
        return self.iid0

    @property
    def iid0(self):
        """A ndarray of the iids. Each iid is a ndarray of two strings (a family ID and a case ID) that identifies an individual.

        :rtype: ndarray (length :attr:`.iid_count`) of ndarray (length 2) of strings

        This property (to the degree practical) reads only iid and sid data from the disk, not SNP value data. Moreover, the iid and sid data is read from file only once.

        :Example:

        >>> from pysnptools.snpreader import Bed
        >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300')
        >>> print snp_on_disk.iid[:3] # print the first three iids
        [['POP1' '0']
         ['POP1' '12']
         ['POP1' '44']]
        """
        return self.row

    @property
    def iid1(self):
        """A ndarray of the iids. Each iid is a ndarray of two strings (a family ID and a case ID) that identifies an individual.

        :rtype: ndarray (length :attr:`.iid_count`) of ndarray (length 2) of strings

        This property (to the degree practical) reads only iid and sid data from the disk, not SNP value data. Moreover, the iid and sid data is read from file only once.

        :Example:

        >>> from pysnptools.snpreader import Bed
        >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300')
        >>> print snp_on_disk.iid[:3] # print the first three iids
        [['POP1' '0']
         ['POP1' '12']
         ['POP1' '44']]
        """
        return self.col

    @property
    def iid_count(self):
        """number of iids

        :rtype: integer

        This property (to the degree practical) reads only iid and sid data from the disk, not SNP value data. Moreover, the iid and sid data is read from file only once.
        """
        assert self.iid0 is self.iid1, "When 'iid_count' is used, iid0 must be the same as iid1"
        return self.iid0_count

    @property
    def iid0_count(self):
        """number of iids

        :rtype: integer

        This property (to the degree practical) reads only iid and sid data from the disk, not SNP value data. Moreover, the iid and sid data is read from file only once.
        """
        return self.row_count

    @property
    def iid1_count(self):
        """number of iids

        :rtype: integer

        This property (to the degree practical) reads only iid and sid data from the disk, not SNP value data. Moreover, the iid and sid data is read from file only once.
        """
        return self.col_count


    @property
    def row_property(self):
        if not hasattr(self,"_row_property"):
            self._row_property = np.empty((self.row_count,0))
        return self._row_property

    @property
    def col_property(self):
        if not hasattr(self,"_col_property"):
            self._col_property = np.empty((self.col_count,0))
        return self._col_property



    #!!check that views always return contiguous memory by default
    def read(self, order='F', dtype=np.float64, force_python_only=False, view_ok=False):
        """Reads the SNP values and returns a :class:`.SnpData` (with :attr:`.SnpData.val` property containing a new ndarray of the SNP values).


        :param order: {'F' (default), 'C', 'A'}, optional -- Specify the order of the ndarray. If order is 'F' (default),
            then the array will be in F-contiguous order (iid-index varies the fastest).
            If order is 'C', then the returned array will be in C-contiguous order (sid-index varies the fastest).
            If order is 'A', then the :attr:`.SnpData.val`
            ndarray may be in any order (either C-, Fortran-contiguous, or even discontiguous).
        :type order: string or None

        :param dtype: {scipy.float64 (default), scipy.float32}, optional -- The data-type for the :attr:`.SnpData.val` ndarray.
        :type dtype: data-type

        :param force_python_only: optional -- If False (default), may use outside library code. If True, requests that the read
            be done without outside library code.
        :type force_python_only: bool


        :param view_ok: optional -- If False (default), allocates new memory for the :attr:`.SnpData.val`'s ndarray. If True,
            if practical and reading from a :class:`SnpData`, will return a new 
            :class:`SnpData` with a ndarray shares memory with the original :class:`SnpData`.
            Typically, you'll also wish to use "order='A'" to increase the chance that sharing will be possible.
            Use these parameters with care because any change to either ndarray (for example, via :meth:`.SnpData.standardize`) will effect
            the others. Also keep in mind that :meth:`read` relies on ndarray's mechanisms to decide whether to actually
            share memory and so it may ignore your suggestion and allocate a new ndarray anyway.
        :type view_ok: bool

        :rtype: :class:`.SnpData`

        Calling the method again causes the SNP values to be re-read and creates a new in-memory :class:`.SnpData` with a new ndarray of SNP values.

        If you request the values for only a subset of the sids or iids, (to the degree practical) only that subset will be read from disk.

        :Example:

        >>> from pysnptools.snpreader import Bed
        >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300') # Specify SNP data on disk
        >>> snpdata1 = snp_on_disk.read() # Read all the SNP data returning a SnpData instance
        >>> print type(snpdata1.val) # The SnpData instance contains a ndarray of the data.
        <type 'numpy.ndarray'>
        >>> subset_snpdata = snp_on_disk[:,::2].read() # From the disk, read SNP values for every other sid
        >>> print subset_snpdata.val[0,0] # Print the first SNP value in the subset
        2.0
        >>> subsub_snpdata = subset_snpdata[:10,:].read(order='A',view_ok=True) # Create an in-memory subset of the subset with SNP values for the first ten iids. Share memory if practical.
        >>> import numpy as np
        >>> # print np.may_share_memory(subset_snpdata.val, subsub_snpdata.val) # Do the two ndarray's share memory? They could. Currently they won't.       
        """
        val = self._read(None, None, order, dtype, force_python_only, view_ok)
        from kerneldata import KernelData
        ret = KernelData(iid0=self.iid0,iid1=self.iid1, val=val, parent_string=str(self))
        return ret

    def iid_to_index(self, list):
        """Takes a list of iids and returns a list of index numbers

        :param list: list of iids
        :type order: list of list of strings

        :rtype: ndarray of int
        
        This method (to the degree practical) reads only iid and sid data from the disk, not SNP value data. Moreover, the iid and sid data is read from file only once.

        :Example:

        >>> from pysnptools.snpreader import Bed
        >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300') # Specify SNP data on disk
        >>> print snp_on_disk.iid_to_index([['POP1','44'],['POP1','12']]) #Find the indexes for two iids.
        [2 1]
        """
        assert self.iid0 is self.iid1, "When 'iid_to_index' is used, iid0 must be the same as iid1"
        return self.iid0_to_index(list)

    def iid0_to_index(self, list):
        """Takes a list of sids and returns a list of index numbers

        :param list: list of sids
        :type list: list of strings

        :rtype: ndarray of int
        
        This method (to the degree practical) reads only iid and sid data from the disk, not SNP value data. Moreover, the iid and sid data is read from file only once.

        :Example:

        >>> from pysnptools.snpreader import Bed
        >>> snp_on_disk = Bed('../tests/datasets/all_chr.maf0.001.N300') # Specify SNP data on disk
        >>> print snp_on_disk.sid_to_index(['1_10','1_13']) #Find the indexes for two sids.
        [2 9]
        """
        return self.row_to_index(list)

    def __getitem__(self, iid_indexer_and_snp_indexer):
        from _subset import _Subset
        try:
            iid0_indexer, iid1_indexer = iid_indexer_and_snp_indexer
        except:
            iid0_indexer = iid_indexer_and_snp_indexer
            iid1_indexer = iid0_indexer

        return _Subset(self, iid0_indexer, iid1_indexer)

    def _assert_iid0_iid1(self):
        assert np.issubdtype(self._row.dtype, str) and len(self._row.shape)==2 and self._row.shape[1]==2, "iid0 should be dtype str, have two dimensions, and the second dimension should be size 2"
        assert np.issubdtype(self._col.dtype, str) and len(self._col.shape)==2 and self._col.shape[1]==2, "iid1 should be dtype str, have two dimensions, and the second dimension should be size 2"


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    #from pysnptools.snpreader import Bed

    #snp_on_disk = Bed('tests/datasets/all_chr.maf0.001.N300') # Specify some data on disk in Bed format
    #subset_snpreader_1 = snp_on_disk[[3,4],:] #index with an array of indexes
    #print subset_snpreader_1.iid_count, subset_snpreader_1.sid_count
    ##2 1015
    #snpdata1 = subset_snpreader_1.read() # read just the two rows of interest from the disk
    #subset_snpreader_2 = snp_on_disk[:,:0:-2] #index with a slice
    #print subset_snpreader_2.iid_count, subset_snpreader_2.sid_count
    ##300 507
    #boolindexes = [s.startswith('23_') for s in snp_on_disk.sid] # create a Boolean index of sids that start '23_'
    #subset_snpreader_3 = snp_on_disk[:,boolindexes] #index with array of Booleans
    #print subset_snpreader_3.iid_count, subset_snpreader_3.sid_count
    ##300 24

    import doctest
    doctest.testmod()
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
    print "done"