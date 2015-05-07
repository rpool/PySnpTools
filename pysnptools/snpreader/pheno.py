import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from snpreader import SnpReader
import pysnptools.util.pheno as pstpheno
import pysnptools.util as pstutil

class Pheno(SnpReader):
    '''
    This is a class that reads into memory from pheno files
    '''

    _ran_once = False

    #!!!this should this have 'missing' as argument
    def __init__(self, input, iid_if_none=None):
        '''
        input    : string of the name of the file or an in-memory dictionary
        '''
        self.input = input #!!!cmk change to _input
        self._iid_if_none = iid_if_none

    def __repr__(self): 
        return "{0}('{1}')".format(self.__class__.__name__,self.input)


    @property
    def row(self):
        self.run_once()
        return self._row

    @property
    def col(self):
        self.run_once()
        return self._col

    @property
    def col_property(self):
        self.run_once()
        return self._col_property

    def run_once(self):
        if self._ran_once:
            return
        self._ran_once = True

        #!!!cmk switch it, so the main code is here rather than in loadPhen
        if isinstance(self.input,str):
            pheno_input = pstpheno.loadPhen(self.input,missing="") #!!!cmkwhat about missing=-9?
        elif self.input is None: #!!!cmk test this
            assert self._iid_if_none is not None, "If input is None then iid_if_none be given"
            pheno_input = {
            'header':np.empty((0),dtype='str'),
            'vals': np.empty((len(self._iid_if_none), 0)),
            'iid': self._iid_if_none
            }

        else:
            pheno_input = self.input

        if len(pheno_input['header']) > 0 and pheno_input['header'][0] is None:
            pheno_input['header'] = ["pheno{0}".format(i) for i in xrange(len(pheno_input['header']))] #!!!cmk move to reader?

        if len(pheno_input['vals'].shape) == 1:
            pheno_input = {
            'header' : pheno_input['header'],
            'vals' : np.reshape(pheno_input['vals'],(-1,1)),
            'iid' : pheno_input['iid']
            }

        self._row = pheno_input['iid']
        self._col = np.array(pheno_input['header'],dtype='str')
        self._col_property = np.empty((len(self._col),3))
        self._col_property.fill(np.nan)
        self._val = pheno_input['vals']

        assert len(self._row) == self._val.shape[0], "Expect # of iids to match number of rows to values"
        assert len(self._col) == self._val.shape[1], "Expect # of sids to match number of columns to values"

        return self

    #def __del__(self):
    #    if self._filepointer != None:  # we need to test this because Python doesn't guarantee that __init__ was fully run
    #        self._filepointer.close()

    def copyinputs(self, copier):
        # doesn't need to self.run_once()
        copier.input(SnpReader.input) #!!!cmk test that this works when input is an inmemory dictionary

    def _read(self, iid_index_or_none, sid_index_or_none, order, dtype, force_python_only, view_ok):
        '''
        Output dictionary:
        'iid' : [N*2] array of family IDs and individual IDs
        'sid' : [S] array rs-numbers or snp identifiers
        'pos' : [S*3] array of positions [chromosome, genetic dist, basepair dist]
        'val' : [N*S] matrix of per iid snp values
        '''
        assert not hasattr(self, 'ind_used'), "A SnpReader should not have a 'ind_used' attribute"
        if order is None:
            order = "F"
        if dtype is None:
            dtype = np.float64
        if force_python_only is None:
            force_python_only = False #!!!cmk this is ignored

        #This could be re-factored to not use so many names
        #!!!cmk similar code appears several places
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


        val = pstutil.sub_matrix(self._val, iid_index_out, sid_index_out, order=order, dtype=dtype)
        
        return val

#!!!cmk add a static write method

if __name__ == "__main__":
    #!!!cmk add test code
    pass