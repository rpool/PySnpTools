import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from pstreader import PstReader
from pstdata import PstData

class _Subset(PstReader):

    def __init__(self, internal, row_indexer, col_indexer):
        '''
        an indexer can be:
             an integer i (same as [i])
             a slice
             a list of integers
             a list of booleans
        '''
        _Subset.static__init__(self, internal, row_indexer, col_indexer)

    @staticmethod
    def static__init__(self, internal, row_indexer, col_indexer):
        self._internal = internal
        self._row_indexer = PstReader._make_sparray_or_slice(row_indexer)
        self._col_indexer = PstReader._make_sparray_or_slice(col_indexer)

    _ran_once = False


    def __repr__(self):
        return _Subset.static__repr__(self) #!!!cmk test this

    def copyinputs(self, copier):
        _Subset.static_copyinput(self,copier)

    @property
    def row(self):
        return _Subset.static_row(self)

    @property
    def col(self):
        return _Subset.static_col(self)

    @property
    def row_property(self):
        return _Subset.static_row_property(self)

    @property
    def col_property(self):
        return _Subset.static_col_property(self)

    def _read(self, row_indexer, col_indexer, order, dtype, force_python_only, view_ok):
        return _Subset.static__read(self, row_indexer, col_indexer, order, dtype, force_python_only, view_ok)

    def run_once(self):
        return _Subset.static_run_once(self)


    @staticmethod
    def static__repr__(self): #!!Should this be __str__ (and elsewhere) because it uses "nice_string" which uses "..." so may be ambiguous?
        s = "{0}[{1},{2}]".format(self._internal,_Subset.static_nice_string(self,self._row_indexer),_Subset.static_nice_string(self,self._col_indexer))
        return s

    _slice_format = {(False,False,False):":",
                     (False,False,True):"::{2}",
                     (False,True,False):":{1}",
                     (False,True,True):":{1}:{2}",
                     (True,False,False):"{0}:",
                     (True,False,True):"{0}::{2}",
                     (True,True,False):"{0}:{1}",
                     (True,True,True):"{0}:{1}:{2}"}

    @staticmethod
    def static_nice_string(self, some_slice):
        if isinstance(some_slice,slice):
            return _Subset._slice_format[(some_slice.start is not None, some_slice.stop is not None, some_slice.step is not None)].format(some_slice.start, some_slice.stop, some_slice.step)
        elif len(some_slice) == 1:
            return str(some_slice[0])
        elif len(some_slice) < 10:
            return "[{0}]".format(",".join([str(i) for i in some_slice]))
        else:
            return "[{0},...]".format(",".join([str(i) for i in some_slice[:10]]))

    @staticmethod
    def static_copyinputs(self, copier):
        self._internal.copyinputs(copier)

    @staticmethod
    def static_row(self):
        self.run_once()
        return self._row

    @staticmethod
    def static_col(self):
        self.run_once()
        return self._col

    @staticmethod
    def static_row_property(self):
        self.run_once()
        return self._row_property

    @staticmethod
    def static_col_property(self):
        self.run_once()
        return self._col_property

    #!!commented out because doesn't guarantee that the shortcut will return with the dtype and order requested.
    #                  Also, didn't handle stacked do-nothing subsets
    #def read(self, order='F', dtype=np.float64, force_python_only=False, view_ok=False):
    #    if view_ok and hasattr(self._internal,"val") and _Subset._is_all_slice(self._row_indexer) and _Subset._is_all_slice(self._col_indexer):
    #        return self._internal
    #    else:
    #        return PstReader.read(self, order, dtype, force_python_only, view_ok)


    # Most _read's support only indexlists or None, but this one supports Slices, too.
    _read_accepts_slices = None
    @staticmethod
    def static__read(self, row_indexer, col_indexer, order, dtype, force_python_only, view_ok):
        self.run_once()

        if hasattr(self._internal,'_read_accepts_slices'):
            composed_row_index_or_none = _Subset.compose_indexer_with_indexer(self._internal.row_count, self._row_indexer, self.row_count, row_indexer)
            composed_col_index_or_none = _Subset.compose_indexer_with_indexer(self._internal.col_count, self._col_indexer, self.col_count, col_indexer)
            val = self._internal._read(composed_row_index_or_none, composed_col_index_or_none, order, dtype, force_python_only, view_ok)
            return val
        else:
            row_index_or_none = PstReader._make_sparray_from_sparray_or_slice(self.row_count, row_indexer)
            composed_row_index_or_none = _Subset.compose_indexer_with_index_or_none(self._internal.row_count, self._row_indexer, self.row_count, row_index_or_none)
            col_index_or_none = PstReader._make_sparray_from_sparray_or_slice(self.col_count, col_indexer)
            composed_col_index_or_none = _Subset.compose_indexer_with_index_or_none(self._internal.col_count, self._col_indexer, self.col_count, col_index_or_none)
            val = self._internal._read(composed_row_index_or_none, composed_col_index_or_none, order, dtype, force_python_only, view_ok)
            return val

    @staticmethod
    def static_run_once(self):
        if self._ran_once:
            return

        self._ran_once = True
        self._row = self._internal.row[self._row_indexer]
        self._col = self._internal.col[self._col_indexer]
        if np.array_equal(self._row,self._col): #When an object is square, keep the row and col the same object.
            self._col = self._row
        self._row_property = self._internal.row_property[self._row_indexer]
        self._col_property = self._internal.col_property[self._col_indexer]
        #self._assert_row_col_pos()

    @staticmethod
    def compose_indexer_with_index_or_none(countA, indexerA, countB, index_or_noneB):
        if _Subset._is_all_slice(indexerA):
            return index_or_noneB

        indexA = PstReader._make_sparray_from_sparray_or_slice(countA, indexerA)

        if _Subset._is_all_slice(index_or_noneB):
            return indexA

        indexAB = indexA[index_or_noneB]

        return indexAB


    @staticmethod
    def compose_indexer_with_indexer(countA, indexerA, countB, indexerB):
        if _Subset._is_all_slice(indexerA):
            return indexerB

        if _Subset._is_all_slice(indexerB):
            return indexerA

        indexA = PstReader._make_sparray_from_sparray_or_slice(countA, indexerA)
        indexB = PstReader._make_sparray_from_sparray_or_slice(countB, indexerB)

        indexAB = indexA[indexB]

        return indexAB
