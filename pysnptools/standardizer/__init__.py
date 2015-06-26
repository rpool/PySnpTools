'''
Define classes such as :class:`Unit` and :class:`Beta` to be used by the :meth:`.SnpData.standardize` method.
'''


import numpy as np
import logging

from pysnptools.standardizer.standardizer import Standardizer
from pysnptools.standardizer.beta import Beta
from pysnptools.standardizer.unit import Unit
from pysnptools.standardizer.identity import Identity
from pysnptools.standardizer.bysidcount import BySidCount
from pysnptools.standardizer.bysqrtsidcount import BySqrtSidCount
from pysnptools.standardizer.diag_K_to_N import DiagKtoN
from pysnptools.standardizer.betatrained import BetaTrained
from pysnptools.standardizer.unittrained import UnitTrained

