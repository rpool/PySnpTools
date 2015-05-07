import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from pstreader import PstReader
import pysnptools.util as pstutil

class PstNpz(PstReader):
    '''
    This is a class that reads into memory from PstNpz files.
    '''

    _ran_once = False

    def __init__(self, pstnpz_filename):
        '''
        filename    : string of the name of the npz file.
        '''
        PstNpz._static__init__(self,pstnpz_filename)

    def __repr__(self): 
        return PstNpz._static__repr__(self)

    @property
    def row(self):
        return PstNpz._static_row(self)

    @property
    def col(self):
        return PstNpz._static_col(self)

    @property
    def row_property(self):
        return PstNpz._static_row_property(self)

    @property
    def col_property(self):
        return PstNpz._static_col_property(self)

    def run_once(self):
        PstNpz._static_run_once(self)

    #def __del__(self):
    #    if self._filepointer != None:  # we need to test this because Python doesn't guarantee that __init__ was fully run
    #        self._filepointer.close()

    def copyinputs(self, copier):
        PstNpz._static_copyinputs(self, copier)

    def _read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok):
        return PstNpz._static_read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok)

    @staticmethod
    def write(data, npz_filename):
        logging.info("Start writing " + npz_filename)
        np.savez(npz_filename, row=data.row, col=data.col, row_property=data.row_property, col_property=data.col_property,val=data.val)
        logging.info("Done writing " + npz_filename)

    @staticmethod
    def _static__init__(self, pstnpz_filename):
        self._npz_filename = pstnpz_filename

    @staticmethod
    def _static__repr__(self): 
        return "{0}('{1}')".format(self.__class__.__name__,self._npz_filename)

    @staticmethod
    def _static_row(self):
        self.run_once()
        return self._row

    @staticmethod
    def _static_col(self):
        self.run_once()
        return self._col

    @staticmethod
    def _static_row_property(self):
        self.run_once()
        return self._row_property

    @staticmethod
    def _static_col_property(self):
        self.run_once()
        return self._col_property


    @staticmethod
    def _static_run_once(self):
        if (self._ran_once):
            return
        self._ran_once = True

        #!!!cmk is this really done without reading 'data'? could mmap support be used?
        with np.load(self._npz_filename) as data: #!! similar code in epistasis
            if len(data.keys()) == 2 and  'arr_0' in data.keys(): #for backwards compatibility
                self._row = data['arr_0']
                self._col = self._row
                self._row_property = np.empty((len(self._row),0))
                self._col_property = np.empty((len(self._col),0))
            else:
                self._row = data['row']
                self._col = data['col']
                if np.array_equal(self._row, self._col): #If it's square, mark it so by making the col and row the same object
                    self._col = self._row
                self._row_property = data['row_property']
                self._col_property = data['col_property']
        #!!!cmk??? self._assert_iid_sid_pos()

        return self

    #def __del__(self):
    #    if self._filepointer != None:  # we need to test this because Python doesn't guarantee that __init__ was fully run
    #        self._filepointer.close()

    @staticmethod
    def _static_copyinputs(self, copier):
        # doesn't need to self.run_once()
        copier.input(self._npz_filename)

    @staticmethod
    def _static_read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok):
        if order is None:
            order = "F"
        if dtype is None:
            dtype = np.float64
        if force_python_only is None:
            force_python_only = False #!!!cmk this is ignore
        #!!!cmk view_ok is also ignored

        #This could be re-factored to not use so many names
        row_count_in = self.row_count
        col_count_in = self.col_count

        if row_index_or_none is not None:
            row_count_out = len(row_index_or_none)
            row_index_out = row_index_or_none
        else:
            row_count_out = row_count_in
            row_index_out = range(row_count_in)

        if col_index_or_none is not None:
            col_count_out = len(col_index_or_none)
            col_index_out = col_index_or_none
        else:
            col_count_out = col_count_in
            col_index_out = range(col_count_in)

        #!!!cmk do we really need to load twice?
        with np.load(self._npz_filename) as data: #!! similar code in epistasis
            if len(data.keys()) == 2 and  'arr_1' in data.keys(): #for backwards compatibility
               val = data['arr_1']
            else:
               val = data['val']

            val = pstutil.sub_matrix(val, row_index_out, col_index_out, order=order, dtype=dtype) #!!!cmk fix so doesn't make a copy of it doesn't need to

        return val



if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    snpreader = Dat(r'../tests/datasets/all_chr.maf0.001.N300.dat')
    snp_matrix = snpreader.read()
    print len(snp_matrix['sid'])
    snp_matrix = snpreader[:,:].read()
    print len(snp_matrix['sid'])
    sid_index_list = snpreader.sid_to_index(['23_9','23_2'])
    snp_matrix = snpreader[:,sid_index_list].read()
    print ",".join(snp_matrix['sid'])
    snp_matrix = snpreader[:,0:10].read()
    print ",".join(snp_matrix['sid'])

    print snpreader.iid_count
    print snpreader.sid_count
    print len(snpreader.pos)

    snpreader2 = snpreader[::-1,4]
    print snpreader.iid_count
    print snpreader2.sid_count
    print len(snpreader2.pos)

    snp_matrix = snpreader2.read()
    print len(snp_matrix['iid'])
    print len(snp_matrix['sid'])

    snp_matrix = snpreader2[5,:].read()
    print len(snp_matrix['iid'])
    print len(snp_matrix['sid'])

    iid_index_list = snpreader2.iid_to_index(snpreader2.iid[::2])
    snp_matrix = snpreader2[iid_index_list,::3].read()
    print len(snp_matrix['iid'])
    print len(snp_matrix['sid'])

    snp_matrix = snpreader[[4,5],:].read()
    print len(snp_matrix['iid'])
    print len(snp_matrix['sid'])

    print snpreader2
    print snpreader[::-1,4]
    print snpreader2[iid_index_list,::3]
    print snpreader[:,sid_index_list]
    print snpreader2[5,:]
    print snpreader[[4,5],:]
