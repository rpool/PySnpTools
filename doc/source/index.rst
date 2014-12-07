################################
:mod:`pysnptools` Documentation
################################

PySnpTools: A library for reading and manipulating genetic data.

:synopsis:

* :mod:`.snpreader`: Efficiently read genetic PLINK formats including \*.bed/bim/fam files. Also, efficiently read *parts* of files and standardize data.

* :mod:`.util`: In one line, intersect and re-order IIDs from :mod:`.snpreader` and other sources. Also, efficiently extract a submatrix from an ndarray.

* :class:`.util.IntRangeSet`: Efficiently manipulate ranges of integers -- for example, genetic position -- with set operators including
  union, intersection, and set difference. 

* :mod:`.util.pheno`: Read the PLINK pheno type file format.

.. automodule:: pysnptools
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:

***********************
:mod:`snpreader` Module
***********************

.. automodule:: pysnptools.snpreader
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:


:class:`snpreader.SnpReader`
+++++++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.SnpReader
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs


:class:`snpreader.Bed`
++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.Bed
    :members:
    :inherited-members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs


:class:`snpreader.snpdata`
++++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.SnpData
    :members:
    :inherited-members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs


:class:`snpreader.ped`
+++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.Ped
    :members:
    :inherited-members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs

:class:`snpreader.dat`
+++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.Dat
    :members:
    :inherited-members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs

:class:`snpreader.hdf5`
+++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.Hdf5
    :members:
    :inherited-members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs


***********************
:mod:`util` Module
***********************

.. automodule:: pysnptools.util
    :members:
    :undoc-members:
	:show-inheritance:

:class:`util.IntRangeSet`
++++++++++++++++++++++++++
.. autoclass:: pysnptools.util.IntRangeSet
  :members:
  :undoc-members:
  :show-inheritance:
  :special-members:
  :exclude-members: __and__, __weakref__,__module__,__dict__, __add__


:mod:`util.pheno`
++++++++++++++++++++++++++
.. automodule:: pysnptools.util.pheno
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
 
****************************
:mod:`standardizer` Module
****************************

.. automodule:: pysnptools.standardizer
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:


:class:`standardizer.Unit`
++++++++++++++++++++++++++
.. autoclass:: pysnptools.standardizer.Unit
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs

:class:`standardizer.Identity`
++++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.standardizer.Identity
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs


:class:`standardizer.Beta`
++++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.standardizer.Beta
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs

:class:`standardizer.BySidCount`
++++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.standardizer.BySidCount
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs

:class:`standardizer.BySqrtSidCount`
++++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.standardizer.BySqrtSidCount
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs

:class:`standardizer.DiagKtoN.py`
++++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.standardizer.DiagKtoN
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs

.. only:: html 

***********************
Indices and Tables
***********************

   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`
