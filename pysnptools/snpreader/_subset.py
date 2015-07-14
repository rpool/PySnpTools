from snpreader import SnpReader
from pysnptools.pstreader._subset import _Subset as PstSubset

class _Subset(PstSubset,SnpReader):
    pass
