import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from pysnptools.snpreader import SnpReader
from pysnptools.snpreader import SnpData
import pysnptools.util.pheno as pstpheno
import pysnptools.util as pstutil
from pysnptools.pstreader import _OneShot

class Pheno(_OneShot, SnpReader):
    '''
    This is a class that reads into memory from pheno files or from in-memory pheno dictionaries
    '''

    def __init__(self, input, iid_if_none=None, missing=None):
        '''
        input    : string of the name of the file or an in-memory dictionary
        '''
        self.filename = input
        self._iid_if_none = iid_if_none
        self.missing = missing

    def _read_pstdata(self):
        #LATER switch it, so the main code is here rather than in loadPhen
        if isinstance(self.filename,str):
            pheno_input = pstpheno.loadPhen(self.filename,missing=self.missing)
        elif self.filename is None:
            assert self._iid_if_none is not None, "If input is None then iid_if_none be given"
            pheno_input = {
            'header':np.empty((0),dtype='str'),
            'vals': np.empty((len(self._iid_if_none), 0)),
            'iid': self._iid_if_none
            }
        else:
            pheno_input = self.filename


        if len(pheno_input['vals'].shape) == 1:
            pheno_input = {
            'header' : pheno_input['header'],
            'vals' : np.reshape(pheno_input['vals'],(-1,1)),
            'iid' : pheno_input['iid']
            }

        if len(pheno_input['header']) > 0 and pheno_input['header'][0] is None:
            pheno_input['header'] = ["pheno{0}".format(i) for i in xrange(len(pheno_input['header']))] #LATER move to reader?
        elif len(pheno_input['header']) == 0:
            pheno_input['header'] = ["pheno{0}".format(i) for i in xrange(pheno_input['vals'].shape[1])]

        row = pheno_input['iid']
        col = np.array(pheno_input['header'],dtype='str')
        col_property = np.empty((len(col),3))
        col_property.fill(np.nan)
        val = pheno_input['vals']

        snpdata = SnpData(iid=row,sid=col,pos=col_property,val=val)
        return snpdata

    @staticmethod
    def write(filename, snpdata, missing='NaN', sep="\t"):
        '''
        !!!cmk
        '''
        with open(filename, 'w') as f:
            for i in xrange(snpdata.iid_count):
                tmpstr = snpdata.iid[i,0] + sep + snpdata.iid[i,1]
                for m in xrange(snpdata.sid_count):
                    v = snpdata.val[i,m]
                    if np.isnan(v):
                        vs = missing
                    else:
                        vs = str(v)
                    tmpstr += sep + vs
                tmpstr += "\n"
                f.write(tmpstr)

