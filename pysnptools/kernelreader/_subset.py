import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from kernelreader import KernelReader
from kerneldata import KernelData
from pysnptools.pstreader._subset import _Subset as PstSubset

class _Subset(KernelReader):

    def __init__(self, internal, row_indexer, col_indexer):
        PstSubset.static__init__(self, internal, row_indexer, col_indexer)

    _ran_once = False

    def __repr__(self):
        return PstSubset.static__repr__(self) #!!!cmk test this

    def copyinputs(self, copier):
        PstSubset.static_copyinput(self,copier)

    @property
    def row(self):
        return PstSubset.static_row(self)

    @property
    def col(self):
        return PstSubset.static_col(self)

    @property
    def row_property(self):
        return PstSubset.static_row_property(self)

    @property
    def col_property(self):
        return PstSubset.static_col_property(self)

    _read_accepts_slices = None
    def _read(self, row_indexer, col_indexer, order, dtype, force_python_only, view_ok):
        return PstSubset.static__read(self, row_indexer, col_indexer, order, dtype, force_python_only, view_ok)

    def run_once(self):
        return PstSubset.static_run_once(self)
