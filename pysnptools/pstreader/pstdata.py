import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from pstreader import PstReader

class PstData(PstReader):
    """  This is a class hold SNP values in-memory along with related iid and sid information.
    It is created by calling the :meth:`.SnpReader.read` method on another :class:`.SnpReader`, for example, :class:`.Bed`.

    See :class:`.SnpReader` for details and examples.
    """
    def __init__(self, row, col, row_property, col_property, val, parent_string="",copyinputs_function=None): #!!!autodoc doesn't generate good doc for this constructor
        #!!!cmk what checks should be done?
        self._row = row #if len(row)>0 else np.array([],dtype=str).reshape(0,2)
        self._col = col #if len(col)>0 else np.array([],dtype=str)
        self._row_property = row_property
        self._col_property = col_property
        #self._pos = pos if len(col)>0 else np.array([],dtype=int).reshape(0,3)

        #self._assert_row_col_pos()

        assert type(val) == np.ndarray, "expect SnpData's val to be a ndarray"
        self.val = val
        self._parent_string = parent_string
        self._std_string_list = []
        if copyinputs_function is not None:
            self.copyinputs = copyinputs_function

    val = None
    """The in-memory data. A numpy.ndarray with dimensions :attr:`.row_count` x :attr:`.col_count`.

    See :class:`.PstReader` for details and examples.
    """

    def __repr__(self):
        if self._parent_string == "":
            if len(self._std_string_list) > 0:
                s = "{0}({1})".format(",".join(self.__class__.__name__,self._std_string_list))
            else:
                s = "{0}()".format(self.__class__.__name__)
        else:
            if len(self._std_string_list) > 0:
                s = "{0}({1},{2})".format(self.__class__.__name__,self._parent_string,",".join(self._std_string_list))
            else:
                s = "{0}({1})".format(self.__class__.__name__,self._parent_string)
        return s

    def copyinputs(self, copier):
        pass

    @property
    def row(self):
        """A ndarray of the iids.

        See :attr:`.SnpReader.iid` for details and examples.
        """
        return self._row

    @property
    def col(self):
        """A ndarray of the sids.

        See :attr:`.SnpReader.sid` for details and examples.
        """
        return self._col

    @property
    def row_property(self):
        """A ndarray of the position information for each sid.

        See :attr:`.SnpReader.pos` for details and examples.
        """
        return self._row_property
    @property
    def col_property(self):
        """A ndarray of the position information for each sid.

        See :attr:`.SnpReader.pos` for details and examples.
        """
        return self._col_property

    #!!Seems like we can't short cut the view_OK this because the .val wouldn't necessarily have the right order and dtype
    #def read(self, order='F', dtype=np.float64, force_python_only=False, view_ok=False):
    #    """creates an in-memory :class:`.SnpData` with a copy of the genotype data
    #    """
    #    if view_ok:
    #        return self
    #    else:
    #        return SnpReader.read(self, order, dtype, force_python_only, view_ok)


    # Most _read's support only indexlists or None, but this one supports Slices, too.
    _read_accepts_slices = None
    def _read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok):
        val, shares_memory = self._apply_sparray_or_slice_to_val(self.val, row_index_or_none, col_index_or_none, order, dtype, force_python_only)
        if shares_memory and not view_ok:
            val = val.copy(order='K')
        return val


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
