import numpy as np
import scipy as sp
import logging
import doctest
import unittest
import os.path
import time

class TestLoader(unittest.TestCase):     

    def test_read(self):
        from pysnptools.pstreader import PstData
        np.random.seed(0)
        row_property=np.array([[1.0,2,2.5],[3,4,4.5],[5,6,6.5]])
        col_property=np.array([[1.0,2,2.5,1],[3,4,4.5,3]])
        pstdata = PstData(row=np.array([[1.0,2],[3,4],[5,6]]), #!!!cmk test 1d and 2d, test various types
                          col=np.array([["A","a"],["B","b"]]),#!!!cmk test 1d and 2d, test various types
                          row_property=row_property,#!!!cmk test 1d and 2d and None
                          col_property=col_property,#!!!cmk test 1d and 2d and None
                          val = np.random.normal(.5,2,size=(3,2)),
                          parent_string="test_read")

        assert pstdata.row_to_index([np.array([3.0,4])])[0] == 1
        assert pstdata.col_to_index([np.array(["A","a"])])[0] == 0
        assert np.array_equal(pstdata[1:,:2].row_property,row_property[1:])
        assert np.array_equal(pstdata[1:,:2].col_property,col_property[:2])


        pstdata2 = pstdata[:2,:2].read()
        np.testing.assert_array_almost_equal(pstdata2.val, pstdata.val[:2,:2], decimal=10)
        pstdata3 = pstdata[[],:].read()
        assert pstdata3.val.shape[0] == 0 and pstdata3.val.shape[1]==2
        pstdata.val = pstdata.val.copy(order='F')
        pstdata2 = pstdata[:2,:2].read()
        np.testing.assert_array_almost_equal(pstdata2.val, pstdata.val[:2,:2], decimal=10)
        pstdata2 = pstdata[:2,:2].read(order='F')
        np.testing.assert_array_almost_equal(pstdata2.val, pstdata.val[:2,:2], decimal=10)
        pstdata2 = pstdata[:2,:2].read(order='A')
        np.testing.assert_array_almost_equal(pstdata2.val, pstdata.val[:2,:2], decimal=10)
        pstdata2 = pstdata[:2,:2].read(force_python_only=True,dtype=None,order='C')
        np.testing.assert_array_almost_equal(pstdata2.val, pstdata.val[:2,:2], decimal=10)
        pstdata2 = pstdata[:2,:2].read(force_python_only=True,dtype='float32',order='C')
        np.testing.assert_array_almost_equal(pstdata2.val, pstdata.val[:2,:2].astype(dtype='float32'), decimal=10)
        pstdata2 = pstdata[:2,:2].read(force_python_only=True,dtype='float32',order=None)
        np.testing.assert_array_almost_equal(pstdata2.val, pstdata.val[:2,:2].astype(dtype='float32'), decimal=10)
        pstdata2 = pstdata[:2,:2].read(force_python_only=True,dtype=None,order='F')
        np.testing.assert_array_almost_equal(pstdata2.val, pstdata.val[:2,:2], decimal=10)
        pstdata4 = pstdata[::,::].read(force_python_only=True)
        np.testing.assert_array_almost_equal(pstdata4.val, pstdata.val, decimal=10)


        logging.info("done with test")


    def test_inputs(self):
        from pysnptools.pstreader import PstData
        np.random.seed(0)
        row_property=np.array([1.0,2,2.5])
        col_property=np.array([1,2,2,1])
        pstdata = PstData(row=np.array([1.0,3,6]), #!!!cmk test 1d and 2d, test various types
                          col=np.array(["Aa","Bb"]),#!!!cmk test 1d and 2d, test various types
                          row_property=row_property,#!!!cmk test 1d and 2d and None
                          col_property=col_property,#!!!cmk test 1d and 2d and None
                          val = np.random.normal(.5,2,size=(3,2)),
                          parent_string="test_read")

        assert pstdata.row_to_index([3])[0] == 1
        assert pstdata.col_to_index(["Aa"])[0] == 0
        assert np.array_equal(pstdata[1:,:2].row_property,row_property[1:])
        assert np.array_equal(pstdata[1:,:2].col_property,col_property[:2])
        logging.info("done with test")


    def test_inputs2(self):
        from pysnptools.pstreader import PstData
        np.random.seed(0)
        row_property=None
        col_property=None
        pstdata = PstData(row=np.array([1.0,3,6]), #!!!cmk test 1d and 2d, test various types
                          col=np.array(["Aa","Bb"]),#!!!cmk test 1d and 2d, test various types
                          row_property=row_property,#!!!cmk test 1d and 2d and None
                          col_property=col_property,#!!!cmk test 1d and 2d and None
                          val = np.random.normal(.5,2,size=(3,2)),
                          parent_string="test_read")

        assert pstdata.row_to_index([3])[0] == 1
        assert pstdata.col_to_index(["Aa"])[0] == 0
        assert np.array_equal(pstdata[1:,:2].row_property,row_property[1:])
        assert np.array_equal(pstdata[1:,:2].col_property,col_property[:2])
        logging.info("done with test")


# We do it this way instead of using doctest.DocTestSuite because doctest.DocTestSuite requires modules to be pickled, which python doesn't allow.
# We need tests to be pickleable so that they can be run on a cluster.
class TestDocStrings(unittest.TestCase):
    pass
    #def test_snpreader(self):
    #    import pysnptools.snpreader.snpreader
    #    old_dir = os.getcwd()
    #    os.chdir(os.path.dirname(os.path.realpath(__file__)))
    #    result = doctest.testmod(pysnptools.snpreader.snpreader)
    #    os.chdir(old_dir)
    #    assert result.failed == 0, "failed doc test: " + __file__

    #def test_bed(self):
    #    import pysnptools.snpreader.bed
    #    old_dir = os.getcwd()
    #    os.chdir(os.path.dirname(os.path.realpath(__file__)))
    #    result = doctest.testmod(pysnptools.snpreader.bed)
    #    os.chdir(old_dir)
    #    assert result.failed == 0, "failed doc test: " + __file__

    #def test_snpdata(self):
    #    import pysnptools.snpreader.snpdata
    #    old_dir = os.getcwd()
    #    os.chdir(os.path.dirname(os.path.realpath(__file__)))
    #    result = doctest.testmod(pysnptools.snpreader.snpdata)
    #    os.chdir(old_dir)
    #    assert result.failed == 0, "failed doc test: " + __file__

    #def test_util(self):
    #    import pysnptools.util
    #    old_dir = os.getcwd()
    #    os.chdir(os.path.dirname(os.path.realpath(__file__))+"/util")
    #    result = doctest.testmod(pysnptools.util)
    #    os.chdir(old_dir)
    #    assert result.failed == 0, "failed doc test: " + __file__


def getTestSuite():
    """
    set up composite test suite
    """
    
    test_suite = unittest.TestSuite([])
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestDocStrings))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestLoader))
    return test_suite

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    suites = getTestSuite()
    r = unittest.TextTestRunner(failfast=True) #!!!cmk
    r.run(suites)
