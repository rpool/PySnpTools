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
from pysnptools.pstreader import PstReader
from pysnptools.snpreader import SnpData
import pysnptools.standardizer as stdizer

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

    #!!!cmk0 Get test_respect_inputs to pass (it has new two-snpreader code
    #!!!cmk0 Think about what happens if there are to snpreaders and they have orderlapping iids
    #!!!cmk0 Think about About if it every makes since to subset a kernel (yes, to get parts. no, subsetting should always be applied to the implict snpreader)
    #!!!cmk0 Think about what intersect_apply should do with: square, xclusive iids rect, overlappyg iid rect
    #!!!cmk0 We want to apply the standaridzation from the training to the test, right? Or we want to standardize test and try together (easy to do because both are right there), right?

    def test_cpp_std(self):

        #Order C vs F
        for order in ['C','F']:
            #32 vs 64
            for dtype in [np.float64,np.float32]:
                #unit vs beta
                for std in [stdizer.Unit(),stdizer.Beta(2,10)]:
                        np.random.seed(0)
                        snp_count = 20
                        snpreader0 = SnpData(iid=[["0","0"],["1","1"],["2","2"]],sid=[str(i) for i in xrange(snp_count)],val=np.array(np.random.randint(3,size=[3,snp_count]),dtype=dtype,order=order))
                        snpreader1 = SnpData(iid=[["3","3"],["4","4"]],sid=[str(i) for i in xrange(snp_count)],val=np.array(np.random.randint(3,size=[2,snp_count]),dtype=dtype,order=order))

                        #has SNC
                        for has_SNC_in_train in [False, True]:
                            if has_SNC_in_train:
                                snpreader0.val[:,1] = 0

                            #missing data
                            for has_missing_data in [False, True]:
                                if has_missing_data:
                                    snpreader0.val[0,2]=np.nan
                                    snpreader1.val[0,2]=np.nan

                                #apply in place or not
                                for apply_in_place in [True,False]:
                                    #gather stats vs not
                                    cppa = snpreader0.read(order=order,dtype=dtype)
                                    pya = snpreader0.read(order=order,dtype=dtype)
                                    stdcppa = std._train_standardizer(cppa,apply_in_place=apply_in_place,force_python_only=False)
                                    stdpya = std._train_standardizer(pya,apply_in_place=apply_in_place,force_python_only=True)
                                    np.testing.assert_array_almost_equal(cppa.val, pya.val, decimal=10 if dtype==np.float64 else 5)

                                    if not apply_in_place:
                                        np.testing.assert_array_almost_equal(cppa.val, snpreader0.val, decimal=10 if dtype==np.float64 else 5)
                                    np.testing.assert_array_almost_equal(stdcppa.stats,stdpya.stats, decimal=10 if dtype==np.float64 else 5)
                                    assert (np.inf in stdcppa.stats[:,1]) == has_SNC_in_train
                                    assert (np.inf in stdpya.stats[:,1]) == has_SNC_in_train

                                    if has_SNC_in_train and apply_in_place:
                                        assert np.array_equal(cppa.val[:,1],np.zeros([cppa.val.shape[0]]))
                                        assert np.array_equal(pya.val[:,1],np.zeros([pya.val.shape[0]]))

                                    if has_missing_data:
                                        if apply_in_place:
                                            assert 0 == cppa.val[0,2]
                                            assert 0 == pya.val[0,2]
                                        else:
                                            assert np.isnan(cppa.val[0,2])
                                            assert np.isnan(pya.val[0,2])
                                        
                                    #uses stats
                                    cppb = snpreader1.read(order=order,dtype=dtype)
                                    pyb = snpreader1.read(order=order,dtype=dtype)
                                    stdcppa._train_standardizer(cppb,apply_in_place=True,force_python_only=False)
                                    stdpya._train_standardizer(pyb,apply_in_place=True,force_python_only=True)
                                    np.testing.assert_array_almost_equal(cppb.val, pyb.val, decimal=10 if dtype==np.float64 else 5)
                                    np.testing.assert_array_almost_equal(stdcppa.stats,stdpya.stats, decimal=10 if dtype==np.float64 else 5) #Make sure we haven't messed up the train stats

                                    if has_SNC_in_train:
                                        assert np.array_equal(cppb.val[:,1],np.zeros([cppb.val.shape[0]]))
                                        assert np.array_equal(pyb.val[:,1],np.zeros([pyb.val.shape[0]]))

                                    if has_missing_data:
                                        assert cppb.val[0,2]==0
                                        assert pyb.val[0,2]==0
        logging.info("done with 'test_cpp_std'")


    def test_demo_kernel(self):

        from sklearn import cross_validation
        from pysnptools.standardizer import Unit

        snps_all = Bed(self.currentFolder + "/../examples/toydata")


        #!!!cmk0 would be nice to get the statistics of standardization for all of bed here and apply to both later

        k_fold = cross_validation.KFold(n=snps_all.iid_count, n_folds=10, random_state=0)
        for train_index, test_index in k_fold:
            snps_train = snps_all[train_index,:]
            snps_test = snps_all[test_index,:]

            # returns a snps_train.iid_count x snps_train.iid_count KernelData
            kernel_train = snps_train.read_kernel(standardizer=Unit(),block_size=1000,force_python_only=False) #force_python_only=True#!!!cmk0
            #     '.val' is the ndarray. 'iid' is the list of iids.
            #     We require that the 'standardizer' be given explicitly because both 'Unit()' and 'Identity()' is reasonable.
            #     By giving a block_size, we tell it read only 10 SNPs at a time, standardizing on the fly.

            # returns a snps_train.iid_count x snps_test.iid_count KernelData
            kernel_test = snps_train.read_kernel(standardizer=Unit(),test=snps_test,block_size=1000,force_python_only=False)#force_python_only=True#!!!Cmk0
            #     By giving a 'reference' we tell it to standardize 'snps_test' according to 'snps_train'.
            #     '.val' is the ndarray. 'iid0' is the list of test iids, 'iid1' is the list of train iids.
            #     By giving a block_size, we tell it read everything from disk again, only 10 SNPs at a time, standardizing on the fly.





        ## Returns a 40 x 40 kernel
        #kernel_0_0 = bed0.read_kernel(standardizer=Unit()) 
        ## '.val' is the ndarray. 'iid' is the list of iids.
        ## We require that the 'standardizer' be given explicitly because both 'Unit()' and 'Identity()' is reasonable. ???????

        ##This will return a 40 x 1 kernel.
        #kernel_0_1 = bed1.read_kernel(reference=bed0,standardizer=Unit())
        ##SHOULD IT BE?
        #kernel_0_1 = bed0.read_kernel(snpreader1=bed1,standardizer=Unit())
        ## OR  1 x 40
        #kernel_1_0 = bed1.read_kernel(reference=bed0,standardizer=Unit())

        ## The standardization of the snps will be based only on the reference.
        ##      In other words, for each snp, the mean and std from the bed0
        ##      will be applied to bed1. bed1's values won't be used to determine the mean or std.
        ## '.val' is the ndarray. 'iid0' is the list of iids from bed0 and 'iid1' will be the iids from bed1.
     




    def test_respect_inputs(self):
        np.random.seed(0)
        for dtype_start,decimal_start in [(np.float32,5),(np.float64,10)]:
            for order_start in ['F','C','A']:
                for snp_count in [20,2]:
                    snpdataX = SnpData(iid=[["0","0"],["1","1"],["2","2"]],sid=[str(i) for i in xrange(snp_count)],val=np.array(np.random.randint(3,size=[3,snp_count]),dtype=dtype_start,order=order_start))
                    for stdx in [stdizer.Beta(1,25),stdizer.Identity(),stdizer.Unit()]:
                        for snpreader0 in [snpdataX,snpdataX[:,1:]]:
                            snpreader1 = snpreader0[1:,:]

                            refdata0 = snpreader0.read()
                            trained_standardizer = refdata0.train_standardizer(apply_in_place=True,standardizer=stdx,force_python_only=True) #!!!cmk00 remove ,force_python_only=True
                            refval0 = refdata0.val.dot(refdata0.val.T)
                            refdata1 = snpreader1.read().standardize(trained_standardizer,force_python_only=True) #!!!cmk00 remove ,force_python_only=True
                            refval1 = refdata0.val.dot(refdata1.val.T)
                            for dtype_goal,decimal_goal in [(np.float32,5),(np.float64,10)]:
                                for order_goal in ['F','C','A']:
                                    k = snpreader0.read_kernel(standardizer=stdx,block_size=1,order=order_goal,dtype=dtype_goal,force_python_only=False) #!!!cmk00 remove ,force_python_only=True)
                                    PstReader._array_properties_are_ok(k.val,order_goal,dtype_goal)
                                    np.testing.assert_array_almost_equal(refval0,k.val, decimal=min(decimal_start,decimal_goal))
                                    k1 = snpreader0.read_kernel(standardizer=stdx,test=snpreader1,block_size=1,order=order_goal,dtype=dtype_goal,force_python_only=False) #!!!cmk00 remove ,force_python_only=True
                                    PstReader._array_properties_are_ok(k1.val,order_goal,dtype_goal)
                                    np.testing.assert_array_almost_equal(refval1,k1.val, decimal=min(decimal_start,decimal_goal))


    def test_respect_nonsquare(self):
        np.random.seed(0)
        snp_count = 20
        snpdata0 = SnpData(iid=[["0","0"],["1","1"],["2","2"]],sid=[str(i) for i in xrange(snp_count)],val=np.array(np.random.randint(3,size=[3,snp_count]),dtype=np.float64))
        snpdata1 = SnpData(iid=[["3","3"],["4","4"]],sid=[str(i) for i in xrange(snp_count)],val=np.array(np.random.randint(3,size=[2,snp_count]),dtype=np.float64))
        #create kernel for snpdata0 x snpdata1
        k01 = snpdata0.read_kernel(standardizer=stdizer.Unit(),test=snpdata1,force_python_only=False) #!!!cmk00 remove ,force_python_only=True)
        k01f = snpdata0.read_kernel(standardizer=stdizer.Unit(),test=snpdata1,order='F',force_python_only=False) #!!!cmk00 remove ,force_python_only=True))
        k0132 = snpdata0.read_kernel(standardizer=stdizer.Unit(),test=snpdata1,dtype=np.float32,force_python_only=False) #!!!cmk00 remove ,force_python_only=True))
        trained_standardizer=snpdata0.train_standardizer(apply_in_place=True,force_python_only=True) #!!!cmk00 remove ,force_python_only=True)
        snpdata1.standardize(trained_standardizer,force_python_only=True) #!!!cmk00 remove ,force_python_only=True))
        refval = snpdata0.val.dot(snpdata1.val.T)
        np.testing.assert_array_almost_equal(refval,k01.val, decimal=10)
        np.testing.assert_array_almost_equal(refval,k01f.val, decimal=10)
        np.testing.assert_array_almost_equal(refval,k0132.val, decimal=5)
        logging.info("done with 'test_respect_nonsquare'")




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
        snpkernel = SnpKernel(snpreader,standardizer=stdizer.Beta())
        s  = str(snpkernel)
        _fortesting_JustCheckExists().input(snpkernel)
        
    def test_npz(self):
        logging.info("in test_npz")
        snpreader = Bed(self.currentFolder + "/../examples/toydata")
        kerneldata1 = snpreader.read_kernel(standardizer=stdizer.Unit())
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
        snpkernel = SnpKernel(snpreader,stdizer.Unit())
        krsub = snpkernel[::2,::2]
        kerneldata1 = krsub.read()
        expected = snpreader[::2,:].kernel(stdizer.Unit())
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
    r = unittest.TextTestRunner(failfast=True) #!!!cmk
    r.run(suites)
