import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from pysnptools.pstreader import PstReader
from pysnptools.pstreader import PstNpz
import pysnptools.util as pstutil
from pysnptools.kernelreader import KernelReader

class KernelNpz(KernelReader):

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
        PstNpz.write(data, npz_filename)

