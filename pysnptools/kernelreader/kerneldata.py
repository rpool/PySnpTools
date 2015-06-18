import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from kernelreader import KernelReader
from pysnptools.pstreader import PstData

class KernelData(KernelReader,PstData):
    """  This is a class hold SNP values in-memory along with related iid and sid information.
    It is created by calling the :meth:`.SnpReader.read` method on another :class:`.SnpReader`, for example, :class:`.Bed`.

    See :class:`.SnpReader` for details and examples.
    """
    def __init__(self, iid=None, iid0=None, iid1=None, val=None, parent_string=""): #!!!autodoc doesn't generate good doc for this constructor
        assert (iid is None) != (iid0 is None and iid1 is None), "Either 'iid' or both 'iid0' 'iid1' must be provided."

        if iid is not None:
            self._row = PstData._fixup_input(iid,empty_creator=lambda ignore:np.empty([0,2],dtype=str))
            self._col = self._row
        else:
            self._row = PstData._fixup_input(iid0,empty_creator=lambda ignore:np.empty([0,2],dtype=str))
            self._col = PstData._fixup_input(iid1,empty_creator=lambda ignore:np.empty([0,2],dtype=str))
        self._row_property = PstData._fixup_input(None,count=len(self._row),empty_creator=lambda count:np.empty([count,0],dtype=str))
        self._col_property = PstData._fixup_input(None,count=len(self._col),empty_creator=lambda count:np.empty([count,0],dtype=str))
        self.val = PstData._fixup_input_val(val,row_count=len(self._row),col_count=len(self._col),empty_creator=lambda row_count,col_count:np.empty([row_count,col_count],dtype=np.float64))

        self._assert_iid0_iid1()
        self._parent_string = parent_string
        self._std_string_list = []

    def standardize(self):
        factor = float(self.iid_count) / np.diag(self.val).sum()
        if abs(factor-1.0)>1e-15:
            self.val *= factor
        return self


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
