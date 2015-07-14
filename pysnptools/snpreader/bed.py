import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from snpreader import SnpReader
from snpdata import SnpData
import math
import warnings

class Bed(SnpReader):
    '''
    A :class:`.SnpReader` for random-access reads of Bed/Bim/Fam files from disk.

    See :class:`.SnpReader` for details and examples.

    The format is described in http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml.

    **Constructor:**
        :Parameters: * **filename** (*string*) -- The \*.bed file to read. The '.bed' suffix is optional. The related \*.bim and \*.fam files will also be read.

    **Methods beyond** :class:`.SnpReader`
    '''
    _ran_once = False
    _filepointer = None

    def __init__(self, filename):
        self.filename = filename

    def __repr__(self): 
        return "{0}('{1}')".format(self.__class__.__name__,self.filename)

    @property
    def row(self):
        """*same as* :attr:`iid`
        """
        self._run_once()
        return self._row

    @property
    def col(self):
        """*same as* :attr:`sid`
        """
        self._run_once()
        return self._col

    @property
    def col_property(self):
        """*same as* :attr:`pos`
        """
        self._run_once()
        return self._col_property

    def _run_once(self):
        if self._ran_once:
            return
        self._ran_once = True

        self._row = SnpReader._read_fam(self.filename,remove_suffix="bed")
        self._col, self._col_property = SnpReader._read_map_or_bim(self.filename,remove_suffix="bed", add_suffix="bim")
        self._assert_iid_sid_pos()

        bedfile = SnpReader._name_of_other_file(self.filename,"bed","bed")
        self._filepointer = open(bedfile, "rb")
        mode = self._filepointer.read(2)
        if mode != 'l\x1b': raise Exception('No valid binary BED file')
        mode = self._filepointer.read(1) #\x01 = SNP major \x00 = individual major
        if mode != '\x01': raise Exception('only SNP-major is implemented')
        logging.info("bed file is open {0}".format(bedfile))

    def __del__(self):
        if self._filepointer != None:  # we need to test this because Python doesn't guarantee that __init__ was fully run
            self._filepointer.close()

    def copyinputs(self, copier):
        # doesn't need to self.run_once() because only uses original inputs
        copier.input(SnpReader._name_of_other_file(self.filename,remove_suffix="bed", add_suffix="bed"))
        copier.input(SnpReader._name_of_other_file(self.filename,remove_suffix="bed", add_suffix="bim"))
        copier.input(SnpReader._name_of_other_file(self.filename,remove_suffix="bed", add_suffix="fam"))


    @staticmethod
    def write(filename, snpdata, force_python_only=False):
        """Writes a :class:`SnpData` to Bed format.

        :param filename: the name of the file to create
        :type filename: string
        :param snpdata: The in-memory data that should be written to disk.
        :type snpdata: :class:`SnpData`

        >>> from pysnptools.snpreader import Pheno, Bed
        >>> import pysnptools.util as pstutil
        >>> snpdata = Pheno('../examples/toydata.phe').read() # Read data from Pheno format
        >>> pstutil.create_directory_if_necessary("tempdir/toydata.bed")
        >>> Bed.write("tempdir/toydata.bed",snpdata)       # Write data in Bed format
        """

        if isinstance(filename,SnpData) and isinstance(snpdata,str): #For backwards compatibility, reverse inputs if necessary
            warnings.warn("write statement should have filename before data to write", DeprecationWarning)
            filename, snpdata = snpdata, filename 


        SnpReader._write_fam(snpdata, filename, remove_suffix="bed")
        SnpReader._write_map_or_bim(snpdata, filename, remove_suffix="bed", add_suffix="bim")

        bedfile = SnpReader._name_of_other_file(filename,remove_suffix="bed", add_suffix="bed")

        if not force_python_only:
            from pysnptools.snpreader import wrap_plink_parser

            if snpdata.val.flags["C_CONTIGUOUS"]:
                order = "C"
            elif snpdata.val.flags["F_CONTIGUOUS"]:
                order = "F"
            else:
                raise Exception("order '{0}' not known, only 'F' and 'C'".format(order))

            if snpdata.val.dtype == np.float64:
                if order=="F":
                    wrap_plink_parser.writePlinkBedFiledoubleFAAA(bedfile, snpdata.iid_count, snpdata.sid_count, snpdata.val)
                else:
                    wrap_plink_parser.writePlinkBedFiledoubleCAAA(bedfile, snpdata.iid_count, snpdata.sid_count, snpdata.val)
            elif snpdata.val.dtype == np.float32:
                if order=="F":
                    wrap_plink_parser.writePlinkBedFilefloatFAAA(bedfile, snpdata.iid_count, snpdata.sid_count, snpdata.val)
                else:
                    wrap_plink_parser.writePlinkBedFilefloatCAAA(bedfile, snpdata.iid_count, snpdata.sid_count, snpdata.val)
            else:
                raise Exception("dtype '{0}' not known, only float64 and float32".format(snpdata.val.dtype))
            
        else:
            with open(bedfile,"wb") as bed_filepointer:
                #see http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml
                bed_filepointer.write(chr(0b01101100)) #magic numbers
                bed_filepointer.write(chr(0b00011011)) #magic numbers
                bed_filepointer.write(chr(0b00000001)) #snp major

                for sid_index in xrange(snpdata.sid_count):
                    if sid_index % 1 == 0:
                        logging.info("Writing snp # {0} to file '{1}'".format(sid_index, filename))

                    col = snpdata.val[:, sid_index]
                    for iid_by_four in xrange(0,snpdata.iid_count,4):
                        vals_for_this_byte = col[iid_by_four:iid_by_four+4]
                        byte = 0b00000000
                        for val_index in xrange(len(vals_for_this_byte)):
                            val = vals_for_this_byte[val_index]
                            if val == 0:
                                code = 0b00
                            elif val == 1:
                                code = 0b10 #backwards on purpose
                            elif val == 2:
                                code = 0b11
                            elif np.isnan(val):
                                code = 0b01 #backwards on purpose
                            else:
                                raise Exception("Can't convert value '{0}' to BED format (only 0,1,2,NAN allowed)".format(val))
                            byte |= (code << (val_index*2))
                        bed_filepointer.write(chr(byte))
        logging.info("Done writing " + filename)


        
    def _read(self, iid_index_or_none, sid_index_or_none, order, dtype, force_python_only, view_ok):
        self._run_once()
        
        if order=='A':
            order='F'

        assert not hasattr(self, 'ind_used'), "A SnpReader should not have a 'ind_used' attribute"

        iid_count_in = self.iid_count
        sid_count_in = self.sid_count

        if iid_index_or_none is not None:
            iid_count_out = len(iid_index_or_none)
            iid_index_out = iid_index_or_none
        else:
            iid_count_out = iid_count_in
            iid_index_out = range(iid_count_in)

        if sid_index_or_none is not None:
            sid_count_out = len(sid_index_or_none)
            sid_index_out = sid_index_or_none
        else:
            sid_count_out = sid_count_in
            sid_index_out = range(sid_count_in)

        if not force_python_only:
            from pysnptools.snpreader import wrap_plink_parser
            val = np.zeros((iid_count_out, sid_count_out), order=order, dtype=dtype)
            bed_fn = SnpReader._name_of_other_file(self.filename,"bed","bed")

            if dtype == np.float64:
                if order=="F":
                    wrap_plink_parser.readPlinkBedFiledoubleFAAA(bed_fn, iid_count_in, sid_count_in, iid_index_out, sid_index_out, val)
                elif order=="C":
                    wrap_plink_parser.readPlinkBedFiledoubleCAAA(bed_fn, iid_count_in, sid_count_in, iid_index_out, sid_index_out, val)
                else:
                    raise Exception("order '{0}' not known, only 'F' and 'C'".format(order));
            elif dtype == np.float32:
                if order=="F":
                    wrap_plink_parser.readPlinkBedFilefloatFAAA(bed_fn, iid_count_in, sid_count_in, iid_index_out, sid_index_out, val)
                elif order=="C":
                    wrap_plink_parser.readPlinkBedFilefloatCAAA(bed_fn, iid_count_in, sid_count_in, iid_index_out, sid_index_out, val)
                else:
                    raise Exception("order '{0}' not known, only 'F' and 'C'".format(order));
            else:
                raise Exception("dtype '{0}' not known, only float64 and float32".format(dtype))
            
        else:
            # An earlier version of this code had a way to read consecutive SNPs of code in one read. May want
            # to add that ability back to the code. 
            # Also, note that reading with python will often result in non-contiguous memory, so the python standardizers will automatically be used, too.       
            logging.warn("using pure python plink parser (might be much slower!!)")
            val = np.zeros(((int(np.ceil(0.25*iid_count_in))*4),sid_count_out),order=order, dtype=dtype) #allocate it a little big
            for SNPsIndex, bimIndex in enumerate(sid_index_out):

                startbit = int(np.ceil(0.25*iid_count_in)*bimIndex+3)
                self._filepointer.seek(startbit)
                nbyte = int(np.ceil(0.25*iid_count_in))
                bytes = np.array(bytearray(self._filepointer.read(nbyte))).reshape((int(np.ceil(0.25*iid_count_in)),1),order='F')

                val[3::4,SNPsIndex:SNPsIndex+1][bytes>=64]=np.nan
                val[3::4,SNPsIndex:SNPsIndex+1][bytes>=128]=1
                val[3::4,SNPsIndex:SNPsIndex+1][bytes>=192]=2
                bytes=np.mod(bytes,64)
                val[2::4,SNPsIndex:SNPsIndex+1][bytes>=16]=np.nan
                val[2::4,SNPsIndex:SNPsIndex+1][bytes>=32]=1
                val[2::4,SNPsIndex:SNPsIndex+1][bytes>=48]=2
                bytes=np.mod(bytes,16)
                val[1::4,SNPsIndex:SNPsIndex+1][bytes>=4]=np.nan
                val[1::4,SNPsIndex:SNPsIndex+1][bytes>=8]=1
                val[1::4,SNPsIndex:SNPsIndex+1][bytes>=12]=2
                bytes=np.mod(bytes,4)
                val[0::4,SNPsIndex:SNPsIndex+1][bytes>=1]=np.nan
                val[0::4,SNPsIndex:SNPsIndex+1][bytes>=2]=1
                val[0::4,SNPsIndex:SNPsIndex+1][bytes>=3]=2
            val = val[iid_index_out,:] #reorder or trim any extra allocation


            #!!LATER this can fail because the trim statement above messes up the order
            #assert(SnpReader._array_properties_are_ok(val, order, dtype)) #!!

        return val


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
