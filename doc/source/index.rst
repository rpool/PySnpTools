################################
:mod:`pysnptools` Documentation
################################

PySnpTools: A library for reading and manipulating genetic data.

:synopsis:

* :mod:`.snpreader`: Efficiently read genetic PLINK formats including \*.bed/bim/fam and phenotype files. Also, efficiently read *parts* of files and standardize data.

* :mod:`.kernelreader`: Efficiently read and manipulate kernel data.

* :mod:`.util`: In one line, intersect and re-order IIDs from :mod:`.snpreader`, :mod:`.kernelreader` and other sources. Also, efficiently extract a submatrix from an ndarray.

* :class:`.util.IntRangeSet`: Efficiently manipulate ranges of integers -- for example, genetic position -- with set operators including
  union, intersection, and set difference. 

* :mod:`.pstreader`: Generalizes :mod:`.snpreader` and :mod:`.kernelreader` (provides the efficiency of numpy arrays with some of the flexibility of pandas)


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

:class:`snpreader.Pheno`
++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.Pheno
    :members:
    :inherited-members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs

:class:`snpreader.SnpData`
++++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.SnpData
    :members:
    :inherited-members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs

:class:`snpreader.Ped`
+++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.Ped
    :members:
    :inherited-members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs

:class:`snpreader.Dat`
+++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.Dat
    :members:
    :inherited-members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs

:class:`snpreader.SnpHdf5`
+++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.SnpHdf5
    :members:
    :inherited-members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs

:class:`snpreader.SnpNpz`
+++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.SnpNpz
    :members:
    :inherited-members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs

***********************
:mod:`kernelreader` Module
***********************

.. automodule:: pysnptools.kernelreader
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:


:class:`kernelreader.KernelReader`
+++++++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelreader.KernelReader
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs


:class:`kernelreader.KernelData`
++++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelreader.KernelData
    :members:
    :inherited-members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs

:class:`kernelreader.Identity`
+++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelreader.Ped
    :members:
    :inherited-members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs


:class:`kernelreader.KernelHdf5`
+++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelreader.KernelHdf5
    :members:
    :inherited-members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs

:class:`kernelreader.KernelNpz`
+++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelreader.KernelNpz
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


***********************
:mod:`pstreader` Module
***********************

.. automodule:: pysnptools.pstreader
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:


:class:`pstreader.PstReader`
+++++++++++++++++++++++++++++
.. autoclass:: pysnptools.pstreader.PstReader
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs


:class:`pstreader.PstData`
++++++++++++++++++++++++++
.. autoclass:: pysnptools.pstreader.PstData
    :members:
    :inherited-members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs

:class:`pstreader.PstHdf5`
+++++++++++++++++++++++++
.. autoclass:: pysnptools.pstreader.PstHdf5
    :members:
    :inherited-members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs

:class:`pstreader.PstNpz`
+++++++++++++++++++++++++
.. autoclass:: pysnptools.pstreader.PstNpz
    :members:
    :inherited-members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs


.. only:: html 

***********************
Indices and Tables
***********************

   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`
