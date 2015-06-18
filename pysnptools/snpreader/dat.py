import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from snpreader import SnpReader
from snpdata import SnpData
import warnings
from pysnptools.pstreader import _OneShot

class Dat(_OneShot,SnpReader):
    '''
    This is a class that reads into memory from DAT/FAM/MAP files.
    '''

    def __init__(self, filename):
        '''
        filename    : string of the name of the Dat file.
        '''
        self.filename = SnpReader._name_of_other_file(filename,remove_suffix="dat", add_suffix="dat")

    def _read_pstdata(self):
        row = SnpReader._read_fam(self.filename,remove_suffix="dat")
        col, col_property = SnpReader._read_map_or_bim(self.filename,remove_suffix="dat", add_suffix="map")
        if len(row)==0 or len(col)==0:
            return SnpData(iid=row,sid=col,pos=col_property,val=np.empty([len(row),len(col)]))
        datfields = pd.read_csv(self.filename,delimiter = '\t',header=None,index_col=False)
        if not np.array_equal(np.array(datfields[0],dtype="string"), col) : raise Exception("Expect snp list in map file to exactly match snp list in dat file")
        del datfields[0]
        del datfields[1]
        del datfields[2]
        assert len(row) == datfields.shape[1], "Expect # iids in fam file to match dat file"
        val = datfields.as_matrix().T
        snpdata = SnpData(iid=row,sid=col,pos=col_property,val=val)
        return snpdata

    def copyinputs(self, copier):
        # doesn't need to self.run_once() because creates name of all files itself
        copier.input(SnpReader._name_of_other_file(self.filename,remove_suffix="dat", add_suffix="dat"))
        copier.input(SnpReader._name_of_other_file(self.filename,remove_suffix="dat", add_suffix="fam"))
        copier.input(SnpReader._name_of_other_file(self.filename,remove_suffix="dat", add_suffix="map"))

    @staticmethod
    def write(filename, snpdata):

        if isinstance(filename,SnpData) and isinstance(snpdata,str): #For backwards compatibility, reverse inputs if necessary
            warnings.warn("write statement should have filename before data to write", DeprecationWarning)
            filename, snpdata = snpdata, filename 

        SnpReader._write_fam(snpdata, filename, remove_suffix="dat")
        SnpReader._write_map_or_bim(snpdata, filename, remove_suffix="dat", add_suffix="map")

        snpsarray = snpdata.val
        with open(filename,"w") as dat_filepointer:
            for sid_index, sid in enumerate(snpdata.sid):
                if sid_index % 1000 == 0:
                    logging.info("Writing snp # {0} to file '{1}'".format(sid_index, filename))
                dat_filepointer.write("{0}\tj\tn\t".format(sid)) #use "j" and "n" as the major and minor allele
                row = snpsarray[:,sid_index]
                dat_filepointer.write("\t".join((str(i) for i in row)) + "\n")
        logging.info("Done writing " + filename)

