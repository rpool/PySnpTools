import numpy as np
import logging
from snpreader import SnpReader
from snpdata import SnpData
import warnings
from pysnptools.pstreader import PstData
from pysnptools.pstreader import _OneShot

class Dense(_OneShot,SnpReader):
    '''
    This is a class that reads into memory from DenseAnsi files. This format looks like:

    var	4006	9570    22899	37676	41236	41978	55159	66828...
    1-10004-rs12354060	22222222...
    1-707348-rs12184279	222222?2...
    1-724325-rs12564807	00000000...
    ...

    where rows are SNPs and columns are observations.
    !!!cmk add example
    '''
    def __init__(self, filename, extract_iid_function=lambda s:("0",s), extract_bim_function=lambda s:("0",s,0,0)):
        '''
        filename    : string of the name of the Dat file.
        '''
        self.filename = filename
        self.extract_iid_function = extract_iid_function
        self.extract_bim_function = extract_bim_function

    def _read_pstdata(self):
        bim_list = []
        val_list_list = []
        with open(self.filename,"r") as fp:
            header = fp.readline()
            iid_string_list = header.strip().split()[1:]
            iid = np.array([self.extract_iid_function(iid_string) for iid_string in iid_string_list],dtype="string")
            val_list = []
            for line_index,line in enumerate(fp):
                if line_index % 1000 == 0:
                    logging.info("reading sid and iid info from line {0} of file '{1}'".format(line_index, self.filename))
                sid_string, rest = line.strip().split()
                assert len(rest) == len(iid)
                bim_list.append(self.extract_bim_function(sid_string))
                val_list = np.array([float(val) if val!="?" else np.NaN for val in rest])
                val_list_list.append(val_list)

        col = np.array([bim[1] for bim in bim_list],dtype='str')
        col_property = np.array([[bim[0],bim[2],bim[3]] for bim in bim_list],dtype=np.float64)

        val = np.zeros((len(iid),len(col)))
        for col_index in xrange(len(col)):
            val[:,col_index] = val_list_list[col_index]

        return PstData(iid,col,val,col_property=col_property,parent_string=self.filename)

    @staticmethod
    def write(filename, snpdata, join_iid_function=lambda iid_pair:iid_pair[1]):

        if isinstance(filename,SnpData) and isinstance(snpdata,str): #For backwards compatibility, reverse inputs if necessary
            warnings.warn("write statement should have filename before data to write", DeprecationWarning)
            filename, snpdata = snpdata, filename 

        snpsarray = snpdata.val
        with open(filename,"w") as filepointer:
            filepointer.write("var"+"\t")
            filepointer.write("\t".join((join_iid_function(iid_pair) for iid_pair in snpdata.iid)) + "\n")

            for sid_index, sid in enumerate(snpdata.sid):
                if sid_index % 1000 == 0:
                    logging.info("Writing snp # {0} to file '{1}'".format(sid_index, filename))
                filepointer.write("{0}\t".format(sid))
                row = snpsarray[:,sid_index]
                filepointer.write("".join((str(int(i)) if i==i else "?" for i in row)) + "\n")
        logging.info("Done writing " + filename)

