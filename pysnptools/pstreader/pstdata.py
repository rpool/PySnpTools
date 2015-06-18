import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from pstreader import PstReader


def _default_empty_creator(count):
    return np.empty([count or 0, 0],dtype=str)

def _default_empty_creator_val(row_count,col_count):
    return np.empty([row_count,col_count],dtype=str)

class PstData(PstReader):
    """  This is a class hold SNP values in-memory along with related iid and sid information.
    It is created by calling the :meth:`.SnpReader.read` method on another :class:`.SnpReader`, for example, :class:`.Bed`.

    See :class:`.SnpReader` for details and examples.
    """
    def __init__(self, row, col, val, row_property=None, col_property=None, parent_string="",copyinputs_function=None): #!!!autodoc doesn't generate good doc for this constructor
        self._row = PstData._fixup_input(row)
        self._col = PstData._fixup_input(col)
        self._row_property = PstData._fixup_input(row_property,count=len(self._row))
        self._col_property = PstData._fixup_input(col_property,count=len(self._col))
        self.val = PstData._fixup_input_val(val,row_count=len(self._row),col_count=len(self._col))
        self._parent_string = parent_string

    val = None
    """The in-memory data. A numpy.ndarray with dimensions :attr:`.row_count` x :attr:`.col_count`.

    See :class:`.PstReader` for details and examples.
    """

    def __eq__(a,b):
        """
        !!!cmk document that equal if have same arrays. The parent string doesn't need to be the same
        """
        try:
            return (np.array_equal(a.row,b.row) and
                    np.array_equal(a.col,b.col) and
                    np.array_equal(a.row_property,b.row_property) and
                    np.array_equal(a.col_property,b.col_property) and
                    np.array_equal(a.val,b.val))
        except:
            return False

    @staticmethod
    def _fixup_input(input,count=None, empty_creator=_default_empty_creator):
        if input is None:
            input = empty_creator(count)
        elif not isinstance(input,np.ndarray):
            input = np.array(input)

        assert count is None or len(input) == count, "Expect length of {0} for input {1}".format(count,input)

        return input

    @staticmethod
    def _fixup_input_val(input,row_count,col_count,empty_creator=_default_empty_creator_val):
        if input is None:
            assert row_count == 0 or col_count == 0, "If val is None, either row_count or col_count must be 0"
            input = _default_empty_creator_val(row_count, col_count)
        elif not isinstance(input,np.ndarray or (input.dtype not in [np.float32,np.float64])):
            input = np.array(input,dtype=np.float64)

        assert len(input.shape)==2 and input.shape[0] == row_count and input.shape[1] == col_count, "Expect val input to have two dimensions of size {0}x{1}".format(row_count,col_count)

        return input



    def __repr__(self):
        if self._parent_string == "":
            return "{0}()".format(self.__class__.__name__)
        else:
            return "{0}({1})".format(self.__class__.__name__,self._parent_string)

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
    _read_accepts_slices = True
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
