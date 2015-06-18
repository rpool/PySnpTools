try:
    import h5py
except:
    pass

import logging
import scipy as sp
from pstreader import PstReader
from pstdata import PstData
import warnings

#!! document the format

class PstHdf5(PstReader):

    _ran_once = False
    h5 = None

    def __init__(self, filename, block_size=5000):
        #!! copy relevant comments from Bed reader
        self.filename=filename
        self.block_size = block_size

    def __repr__(self): 
        return "{0}('{1}')".format(self.__class__.__name__,self.filename) #!!LATER print non-default values, too

    def copyinputs(self, copier):
        copier.input(self.filename)

    @property
    def row(self):
        self.run_once()
        return self._row

    @property
    def col(self):
        self.run_once()
        return self._col

    @property
    def row_property(self):
        self.run_once()
        return self._row_property

    @property
    def col_property(self):
        self.run_once()
        return self._col_property

    def _find_vocab(self):
        vocab_list = [['row','col','val','row_property','col_property'],['iid','sid','val',None,'pos'],['iid','rs','snps',None,'pos']]

        for vocab in vocab_list:
            if all((key is None or key in self.h5) for key in vocab):
                return vocab
        raise Exception("Don't know how to read HDF5 with these keys: " + ",".join(self.h5.iterkeys()))


    def run_once(self):
        if self._ran_once:
            return
        try:
            self.h5 = h5py.File(self.filename, "r")
        except IOError, e:
            raise IOError("Missing or unopenable file '{0}' -- Native error message: {1}".format(self.filename,e))

        row_key,col_key,val_key,row_property_key,col_property_key = self._find_vocab()

        self._row = PstData._fixup_input(self.h5[row_key])
        self._col = PstData._fixup_input(self.h5[col_key])
        self._row_property = PstData._fixup_input(self.h5[row_property_key] if row_property_key else None,count=len(self._row))  #Extra "if ... else" for backwards compatibility.
        self._col_property = PstData._fixup_input(self.h5[col_property_key],count=len(self._col))
        self.val_in_file = self.h5[val_key]

        self.is_col_major = None
        if "col-major" in self.val_in_file.attrs:
            self.is_col_major = self.val_in_file.attrs["col-major"]
        elif "SNP-major" in self.val_in_file.attrs:
            self.is_col_major = self.val_in_file.attrs["SNP-major"]
        assert self.is_col_major is not None, "In Hdf5 the 'val' matrix must have a Boolean 'col-major' (or 'SNP-major') attribute"

        S_original = len(self._col)
        N_original = len(self._row)
        if self.is_col_major:
            if not self.val_in_file.shape == (S_original, N_original) : raise Exception("In Hdf5, the val matrix dimensions don't match those of 'row' and 'col'")
        else:
            if not self.val_in_file.shape == (N_original, S_original) : raise Exception("In Hdf5, the val matrix dimensions don't match those of 'row' and 'col'")

        self._ran_once = True


    @staticmethod
    def _is_sorted_without_repeats(list):
        if len(list) < 2:
            return True
        for i in xrange(1,len(list)):
            if not list[i-1] < list[i]:
                return False
        return True
    

    def __del__(self):
        if self.h5 != None:  # we need to test this because Python doesn't guarantee that __init__ was fully run
            self.h5.close()

    def read_direct(self, val, selection=sp.s_[:,:]):
        if self.is_col_major:
            selection = tuple(reversed(selection))

        if val.flags["F_CONTIGUOUS"]:
            self.val_in_file.read_direct(val.T,selection)
        else:
            self.val_in_file.read_direct(val,selection)

    def create_block(self, block_size, order, dtype):
        matches_order = self.is_col_major == (order =="F")
        opposite_order = "C" if order == "F" else "F"
        if matches_order:
            return sp.empty([len(self._row),block_size], dtype=dtype, order=order)
        else:
            return sp.empty([len(self._row),block_size], dtype=dtype, order=opposite_order)

    def _read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok):
        self.run_once()

        opposite_order = "C" if order == "F" else "F"

        if row_index_or_none is not None:
            row_index_count = len(row_index_or_none)
            row_index_list = row_index_or_none
            row_is_sorted = PstHdf5._is_sorted_without_repeats(row_index_list)
        else:
            row_index_count = self.row_count
            row_index_list = range(self.row_count)
            row_is_sorted = True

        if col_index_or_none is not None:
            col_index_count = len(col_index_or_none)
            col_index_list = col_index_or_none
        else:
            col_index_count = self.col_count
            col_index_list = range(self.col_count)
        #Check if snps and iids indexes are in order and in range
        col_are_sorted = PstHdf5._is_sorted_without_repeats(col_index_list)

        val = sp.empty([row_index_count, col_index_count], dtype=dtype, order=order)

        matches_order = self.is_col_major == (order=="F")
        is_simple = not force_python_only and row_is_sorted and col_are_sorted and matches_order #If 'is_simple' may be able to use a faster reader

        # case 0 -- zero elements in val
        if row_index_count == 0 or col_index_count == 0:
            pass

        # case 1 - all snps & all ids requested
        elif is_simple and col_index_count == self.col_count and row_index_count == self.row_count:
            self.read_direct(val)

        # case 2 - some snps and all ids
        elif is_simple and row_index_count == self.row_count:
            self.read_direct(val, sp.s_[:,col_index_list])

        # case 3 all snps and some ids
        elif is_simple and col_index_count == self.col_count:
            self.read_direct(val, sp.s_[row_index_list,:])

        # case 4 some snps and some ids -- use blocks
        else:
            block_size = min(self.block_size, col_index_count)
            block = self.create_block(block_size, order, dtype)

            if not col_are_sorted:
                col_index_index_list = sp.argsort(col_index_list)
                col_index_list_sorted = col_index_list[col_index_index_list]
            else:
                col_index_index_list = sp.arange(col_index_count)
                col_index_list_sorted = col_index_list

            for start in xrange(0, col_index_count, block_size):
                #print start
                end = min(start+block_size,col_index_count)
                if end-start < block_size:  #On the last loop, the buffer might be too big, so make it smaller
                    block = self.create_block(end-start, order, dtype)
                col_index_list_forblock = col_index_list_sorted[start:end]
                col_index_index_list_forblock = col_index_index_list[start:end]
                self.read_direct(block, sp.s_[:,col_index_list_forblock])
                val[:,col_index_index_list_forblock] = block[row_index_list,:]

        #!!LATER does this test work when the size is 1 x 1 and order if F? iid_index_or_none=[0], sid_index_or_none=[1000] (based on test_blocking_hdf5)
        has_right_order = (order=="C" and val.flags["C_CONTIGUOUS"]) or (order=="F" and val.flags["F_CONTIGUOUS"])
        assert val.shape == (row_index_count, col_index_count) and val.dtype == dtype and has_right_order
        return val




    @staticmethod
    def write(filename, pstdata, dtype='f8',col_major=True):

        if isinstance(filename,PstData) and isinstance(pstdata,str): #For backwards compatibility, reverse inputs if necessary
            warnings.warn("write statement should have filename before data to write", DeprecationWarning)
            filename, pstdata = pstdata, filename 

        if not isinstance(dtype, str) or len(dtype) != 2 or dtype[0] != 'f' : raise Exception("Expect dtype to start with 'f', e.g. 'f4' for single, 'f8' for double")
        val = (pstdata.val.T) if col_major else pstdata.val

        with h5py.File(filename, "w") as h5:
            h5.create_dataset('row', data=pstdata.row)
            h5.create_dataset('col', data=pstdata.col)
            h5.create_dataset('row_property', data=pstdata.row_property)
            h5.create_dataset('col_property', data=pstdata.col_property)
            h5.create_dataset('val', data=val,dtype=dtype,shuffle=True)#compression="gzip", doesn't seem to work with Anaconda
            h5['val'].attrs["col-major"] = col_major

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    #snpreader = Hdf5(r'../tests/datasets/all_chr.maf0.001.N300.hdf5')
    #snp_matrix = snpreader.read()
    #print len(snp_matrix['sid'])
    #snp_matrix = snpreader[:,:].read()
    #print len(snp_matrix['sid'])
    #sid_index_list = snpreader.sid_to_index(['23_9','23_2'])
    #snp_matrix = snpreader[:,sid_index_list].read()
    #print ",".join(snp_matrix['sid'])
    #snp_matrix = snpreader[:,0:10].read()
    #print ",".join(snp_matrix['sid'])

    #print snpreader.iid_count
    #print snpreader.sid_count
    #print len(snpreader.pos)

    #snpreader2 = snpreader[::-1,4]
    #print snpreader.iid_count
    #print snpreader2.sid_count
    #print len(snpreader2.pos)

    #snp_matrix = snpreader2.read()
    #print len(snp_matrix['iid'])
    #print len(snp_matrix['sid'])

    #snp_matrix = snpreader2[5,:].read()
    #print len(snp_matrix['iid'])
    #print len(snp_matrix['sid'])

    #iid_index_list = snpreader2.iid_to_index(snpreader2.iid[::2])
    #snp_matrix = snpreader2[iid_index_list,::3].read()
    #print len(snp_matrix['iid'])
    #print len(snp_matrix['sid'])

    #snp_matrix = snpreader[[4,5],:].read()
    #print len(snp_matrix['iid'])
    #print len(snp_matrix['sid'])

    #print snpreader2
    #print snpreader[::-1,4]
    #print snpreader2[iid_index_list,::3]
    #print snpreader[:,sid_index_list]
    #print snpreader2[5,:]
    #print snpreader[[4,5],:]


    #import doctest
    #doctest.testmod()
