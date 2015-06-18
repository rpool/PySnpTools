import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from snpreader import SnpReader
from snpdata import SnpData
import numpy as np
import warnings
from pysnptools.pstreader import _OneShot

class Ped(_OneShot,SnpReader):
    '''
    This is a class that does a Ped file. For examples of its use see its 'read' method. #!!LATER update comments
    #!!!cmk mention that 012 may become 210
    '''

    def __init__(self, filename, missing = '0'):
        '''
            filename    : string of the filename of the ped file
            missing         : string indicating a missing genotype (default '0')
        '''

        self.filename = SnpReader._name_of_other_file(filename,remove_suffix="ped", add_suffix="ped")
        self.missing = missing

    def _read_pstdata(self):
        col, col_property = SnpReader._read_map_or_bim(self.filename,remove_suffix="ped", add_suffix="map")
        ped = np.loadtxt(self.filename, dtype='str', comments=None)
        row = ped[:,0:2]
        snpsstr = ped[:,6::]
        inan=snpsstr==self.missing
        snps = np.zeros((snpsstr.shape[0],snpsstr.shape[1]/2))
        for i in xrange(snpsstr.shape[1]/2):
            snps[inan[:,2*i],i]=np.nan
            vals=snpsstr[~inan[:,2*i],2*i:2*(i+1)]
            snps[~inan[:,2*i],i]+=(vals==vals[0,0]).sum(1)
        snpdata = SnpData(iid=row,sid=col,pos=col_property,val=snps)
        return snpdata

    def copyinputs(self, copier):
        # doesn't need to self.run_once() because only uses original inputs
        copier.input(self.filename)
        copier.input(SnpReader._name_of_other_file(self.filename,remove_suffix="ped", add_suffix="map"))


    @staticmethod
    def write(filename, snpdata):

        if isinstance(filename,SnpData) and isinstance(snpdata,str): #For backwards compatibility, reverse inputs if necessary
            warnings.warn("write statement should have filename before data to write", DeprecationWarning)
            filename, snpdata = snpdata, filename 

        SnpReader._write_map_or_bim(snpdata, filename, remove_suffix="ped", add_suffix="map")

        # The PED file is a white-space (space or tab) delimited file: the first six columns are mandatory:
        # Family ID
        # Individual ID
        # Paternal ID
        # Maternal ID
        # Sex (1=male; 2=female; other=unknown)
        # Phenotype

        pedfile = SnpReader._name_of_other_file(filename, remove_suffix="ped", add_suffix="ped")
        with open(pedfile,"w") as ped_filepointer:
            for iid_index, iid_row in enumerate(snpdata.iid):
                ped_filepointer.write("{0} {1} 0 0 0 0".format(iid_row[0],iid_row[1]))
                row = snpdata.val[iid_index,:]
                for sid_index, val in enumerate(row):
                    if val == 0:
                        s = "A A"
                    elif val == 1:
                        s = "A G"
                    elif val == 2:
                        s = "G G"
                    elif np.isnan(val):
                        s = "0 0"
                    else:
                        raise Exception("Expect values for ped file to be 0,1,2, or NAN. Instead, saw '{0}'".format(val))
                    ped_filepointer.write("\t"+s)
                ped_filepointer.write("\n")
