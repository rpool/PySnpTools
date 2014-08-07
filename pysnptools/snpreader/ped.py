import numpy as SP
import subprocess, sys, os.path
from itertools import *
from pysnptools.pysnptools.altset_list import *
import pandas as pd
import logging
from snpreader import SnpReader
import numpy as np

#!!!cmk05192014 LATER: Make Ped and/or Dat reader not be all in memory
class Ped(SnpReader):
    '''
    This is a class that does a Ped file. For examples of its use see its 'read' method. #!!!cmk05192014 LATER update comments
    '''

    _ran_once = False
    _filepointer = None

    def __init__(self, basefilename, missing = '0'):
        '''
            basefilename    : string of the basename of [basename].ped and [basename].map
            missing         : string indicating a missing genotype (default '0')
        '''

        self.basefilename = basefilename
        self.missing = missing

    def __repr__(self): 
        missing_string = "" if self.missing == '0' else ",missing='{0}'".format(self.missing)
        return "{0}('{1}'{2})".format(self.__class__.__name__,self.basefilename,missing_string)

    @property
    def iid(self):
        self.run_once()
        return self._iid

    @property
    def sid(self):
        self.run_once()
        return self._sid

    @property
    def pos(self):
        self.run_once()
        return self._pos

    def run_once(self):
        if (self._ran_once):
            return
        self._ran_once = True

        pedfile = self.basefilename+".ped"
        mapfile = self.basefilename+".map"
        map = SP.loadtxt(mapfile,dtype = 'str',comments=None)

        self._sid = map[:,1]
        self._pos = SP.array(map[:,(0,2,3)],dtype = 'float')
        map = None

        ped = SP.loadtxt(pedfile,dtype = 'str',comments=None)
        self._iid = ped[:,0:2]
        snpsstr = ped[:,6::]
        inan=snpsstr==self.missing
        self._snps = SP.zeros((snpsstr.shape[0],snpsstr.shape[1]/2))
        for i in xrange(snpsstr.shape[1]/2):
            self._snps[inan[:,2*i],i]=SP.nan
            vals=snpsstr[~inan[:,2*i],2*i:2*(i+1)]
            self._snps[~inan[:,2*i],i]+=(vals==vals[0,0]).sum(1)
    def copyinputs(self, copier):
        # doesn't need to self.run_once() because only uses original inputs
        copier.input(self.basefilename + ".ped")
        copier.input(self.basefilename + ".map")


    def _read(self, iid_index_or_none, sid_index_or_none, order, dtype, force_python_only, view_ok):
        '''
        Input: a snp_set. Choices include

        Output dictionary: #!!!cmk05192014 LATER update this
        'iid' : [N*2] array of family IDs and individual IDs
        'sid' : [S] array rs-numbers or snp identifiers
        'pos' : [S*3] array of positions [chromosome, genetic dist, basepair dist]
        'val' : [N*S] matrix of per iid snp values

        Examples:

        >>> bed = Bed(r'../../demo/all_chr.maf0.001.N300')
        ... ret = bed.read()
        ... len(ret['sid'])
        ... ret = bed.read(AllSnps())
        ... len(ret['sid'])
        ... ret = bed.read(SnpAndSetName('someset',['23_9','23_2']))
        ... ",".join(ret['sid'])
        ... ret = bed.read(PositionRange(0,10))
        ... ",".join(ret['sid'])
        Loading fam file ../../demo/all_chr.maf0.001.N300.fam
        Loading bim file ../../demo/all_chr.maf0.001.N300.bim
        bed file is open ../../demo/all_chr.maf0.001.N300.bed
        1015
        1015
        '23_9,23_2'
        '1_12,1_34,1_10,1_35,1_28,1_25,1_36,1_39,1_4,1_13'
        closing bed file


        >>> altset_list1 = SnpAndSetNameCollection(r'../../demo/set_input.small.txt') # get the list of snpsets defined in the file
        Reading ../../demo/set_input.small.txt
        >>> altset_list2 = Subset(altset_list1,['set1','set5'])                       # only use a subset of those snpsets
        >>> altset_list3 = MinMaxSetSize(altset_list2, minsetsize=2, maxsetsize=15)   # only use the subset of subsets that contain between 2 & 15 snps (inclusive)
        >>> bed = Bed(r'../../demo/all_chr.maf0.001.N300')
        ... altsetlist_plusbed = altset_list3.addbed(bed)                        # apply altset_list3 to this bed file
        ... len(altsetlist_plusbed)                                              # tell how many snpsets there will be
        ... for snpset_plusbed in altsetlist_plusbed:
        ...      str(snpset_plusbed)                                             # the name of the snpset
        ...      len(snpset_plusbed)                                             # the number of snps in the snpset
        ...      ret = snpset_plusbed.read()
        ...      ",".join(ret['sid'])
        Loading fam file ../../demo/all_chr.maf0.001.N300.fam
        Loading bim file ../../demo/all_chr.maf0.001.N300.bim
        bed file is open ../../demo/all_chr.maf0.001.N300.bed
        1
        'set5'
        13
        '5_12,5_28,5_32,5_5,5_11,5_1,5_9,5_3,5_19,5_7,5_21,5_15,5_23'
        closing bed file

        '''
        self.run_once()
        assert not hasattr(self, 'ind_used'), "A SnpReader should not have a 'ind_used' attribute"

        if order is None:
            order = "F"
        if dtype is None:
            dtype = SP.float64
        if force_python_only is None:
            force_python_only = False

        val, shares_memory = self._apply_sparray_or_slice_to_val(self._snps, iid_index_or_none, sid_index_or_none, order, dtype, force_python_only)
        if shares_memory:
            val = val.copy(order='K')
        return val

    @staticmethod
    def write(snpdata, basefilename):

        pedfile = basefilename + ".ped"
        mapfile = basefilename + ".map"

        #!!!cmk05192014 LATER this could be made faster

        with open(mapfile,"w") as map_filepointer:
            for sid_index, sid in enumerate(snpdata.sid):
                posrow = snpdata.pos[sid_index]
                map_filepointer.write("{0}\t{1}\t{2}\t{3}\n".format(posrow[0], sid, posrow[1], posrow[2]))

         # The PED file is a white-space (space or tab) delimited file: the first six columns are mandatory:
         # Family ID
         # Individual ID
         # Paternal ID
         # Maternal ID
         # Sex (1=male; 2=female; other=unknown)
         # Phenotype

        with open(pedfile,"w") as ped_filepointer:
            for iid_index, iid_row in enumerate(snpdata.iid):
                ped_filepointer.write("{0} {1} 0 0 0 0".format(iid_row[0],iid_row[1]))
                col = snpdata.val[iid_index,:]
                for sid_index, val in enumerate(col):
                    if val == 0:
                        s = "A A"
                    elif val == 1:
                        s = "A G"
                    elif val == 2:
                        s = "G G"
                    elif val == SP.nan:
                        s = "0 0"
                    else:
                        raise Exception("Expect values for ped file to be 0,1,2, or NAN. Instead, saw '{0}'".format(val))
                    ped_filepointer.write("\t"+s)
                ped_filepointer.write("\n")

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    snpreader = Ped(r'../../tests/datasets/all_chr.maf0.001.N300')
    snp_matrix = snpreader.read()

    #from snpreader.hdf5 import Hdf5
    #Hdf5.write(snp_matrix, r'../../tests/datasets/all_chr.maf0.001.N300.hdf5')

    #from snpreader.dat import Dat
    #Dat.write(snp_matrix, r'../../tests/datasets/all_chr.maf0.001.N300.dat')


    print len(snp_matrix['sid'])
    snp_matrix = snpreader[:,:].read()
    print len(snp_matrix['sid'])
    sid_index_list = snpreader.sid_to_index(['23_9','23_2'])
    snp_matrix = snpreader[:,sid_index_list].read()
    print ",".join(snp_matrix['sid'])
    snp_matrix = snpreader[:,0:10].read()
    print ",".join(snp_matrix['sid'])

    print snpreader.iid_count
    print snpreader.sid_count
    print len(snpreader.pos)

    snpreader2 = snpreader[::-1,4]
    print snpreader.iid_count
    print snpreader2.sid_count
    print len(snpreader2.pos)

    snp_matrix = snpreader2.read()
    print len(snp_matrix['iid'])
    print len(snp_matrix['sid'])

    snp_matrix = snpreader2[5,:].read()
    print len(snp_matrix['iid'])
    print len(snp_matrix['sid'])

    iid_index_list = snpreader2.iid_to_index(snpreader2.iid[::2])
    snp_matrix = snpreader2[iid_index_list,::3].read()
    print len(snp_matrix['iid'])
    print len(snp_matrix['sid'])

    snp_matrix = snpreader[[4,5],:].read()
    print len(snp_matrix['iid'])
    print len(snp_matrix['sid'])

    print snpreader2
    print snpreader[::-1,4]
    print snpreader2[iid_index_list,::3]
    print snpreader[:,sid_index_list]
    print snpreader2[5,:]
    print snpreader[[4,5],:]


    #!!!cmk05192014 LATER
    #import doctest
    #doctest.testmod()
