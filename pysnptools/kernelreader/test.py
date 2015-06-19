import numpy as np
import scipy as sp
import logging
import doctest
import unittest
import os.path
import time

from pysnptools.kernelreader import *
from pysnptools.snpreader import Bed
from pysnptools.util import create_directory_if_necessary
import pysnptools.standardizer as std

class _fortesting_JustCheckExists(object): #Implements ICopier

    def __init__(self,doPrintOutputNames=False):
        self.doPrintOutputNames = doPrintOutputNames
    
    def input(self,item):
        if isinstance(item, str):
            if not os.path.exists(item): raise Exception("Missing input file '{0}'".format(item))
        elif hasattr(item,"copyinputs"):
            item.copyinputs(self)
        # else -- do nothing

    def output(self,item):
        if isinstance(item, str):
            if not os.path.exists(item): raise Exception("Missing output file '{0}'".format(item))
            if self.doPrintOutputNames:
                print item
        elif hasattr(item,"copyoutputs"):
            item.copyoutputs(self)
        # else -- do nothing




class TestLoader(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.currentFolder = os.path.dirname(os.path.realpath(__file__))

    def test_fail(self):
        did_fail = True
        try:
            kd = KernelData(iid=[["0"],["1"],["2"]],val=[[1,2,3],[4,5,6],[7,8,9]]) #Wrong iid shape
            did_fail = False
        except:
            pass
        assert did_fail, "The constructor should fail because the iid is the wrong shape"

    def test_kernel2(self):
        logging.info("in kernel2")
        kd = KernelData(iid=[["0","0"],["1","1"],["2","2"]],val=[[1,2,3],[4,5,6],[7,8,9]])
        assert np.array_equal(kd.iid_to_index([["1","1"],["2","2"]]),np.array([1,2]))
        assert np.array_equal(kd.iid0_to_index([["1","1"],["2","2"]]),np.array([1,2]))
        kd = kd.standardize()
        assert np.abs(np.diag(kd.val).sum()-3)<1e-7
        assert kd.iid1_count == 3

    def test_snp_kernel2(self):
        logging.info("in test_snp_kernel2")
        snpreader = Bed(self.currentFolder + "/../examples/toydata")
        snpkernel = SnpKernel(snpreader,standardizer=std.Beta())
        s  = str(snpkernel)
        _fortesting_JustCheckExists().input(snpkernel)
        
    def test_npz(self):
        logging.info("in test_npz")
        snpreader = Bed(self.currentFolder + "/../examples/toydata")
        kerneldata1 = snpreader.read_kernel(std.Unit())
        s = str(kerneldata1)
        output = "tempdir/kernelreader/toydata.kernel.npz"
        create_directory_if_necessary(output)
        KernelNpz.write(output,kerneldata1)
        kernelreader2 = KernelNpz(output)
        kerneldata2 = kernelreader2.read()
        np.testing.assert_array_almost_equal(kerneldata1.val, kerneldata2.val, decimal=10)
        logging.info("done with test")

    def test_subset(self):
        logging.info("in test_subset")
        snpreader = Bed(self.currentFolder + "/../examples/toydata")
        snpkernel = SnpKernel(snpreader,std.Unit())
        krsub = snpkernel[::2,::2]
        kerneldata1 = krsub.read()
        expected = snpreader[::2,:].kernel(std.Unit())
        np.testing.assert_array_almost_equal(kerneldata1.val, expected, decimal=10)

        krsub2 = snpkernel[::2]
        kerneldata2 = krsub2.read()
        np.testing.assert_array_almost_equal(kerneldata2.val, expected, decimal=10)
        logging.info("done with test")

    def test_identity(self):
        logging.info("in test_identity")
        snpreader = Bed(self.currentFolder + "/../examples/toydata")
        assert snpreader.iid is snpreader.row
        kid = Identity(snpreader.row)
        assert np.array_equal(kid.row,kid.iid) and np.array_equal(kid.iid,kid.iid0) and np.array_equal(kid.iid0,kid.iid1) and np.array_equal(kid.iid1, kid.col)
        np.testing.assert_array_almost_equal(kid.read().val, np.identity(snpreader.iid_count))

        #subset
        sub1 = kid[1:5]
        np.testing.assert_array_almost_equal(sub1.read().val, np.identity(4))

        #zero size
        sub2 = kid[0:0]
        np.testing.assert_array_almost_equal(sub2.read().val, np.identity(0))

    def test_identity_sub(self):
        logging.info("in test_identity_sub")
        snpreader = Bed(self.currentFolder + "/../examples/toydata")
        assert snpreader.iid is snpreader.row
        kid = Identity(snpreader.row)
        sub3 = kid[::2,1:5]
        expected = np.identity(kid.iid_count)[::2,:][:,1:5]
        np.testing.assert_array_almost_equal(sub3.read().val,expected)


        logging.info("done with test")


def getTestSuite():
    """
    set up composite test suite
    """
    
    test_suite = unittest.TestSuite([])
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestLoader))
    return test_suite

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    suites = getTestSuite()
    r = unittest.TextTestRunner(failfast=False)
    r.run(suites)
