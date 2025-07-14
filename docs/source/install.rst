Installation
=================

Requirements
------------
Before installing Magnon, please make sure your Python version and packages meet the requirements for

* Python (version 3.13 or newer)
* `NumPy <https://numpy.org/doc/stable/index.html>`_ - for vectorised array handling (version 2.1 or newer)
* `Matplotlib <https://matplotlib.org/stable/>`_ - for plotting (version 3.9.2 or newer)
* `Atomic Simulation Environment <https://wiki.fysik.dtu.dk/ase/index.html>`_ - for handling crystal structure data (version 3.23 or newer)
* `Spglib <https://spglib.readthedocs.io/en/latest/>`_ - to get crystal symmetries (version 2.5 or newer)

all of which are available through the standard package repositories (such as PyPI).

PyPI (pip install)
------------------

Magnon is most easily installed using the Python Package Index using:

.. code-block:: bash

   pip install --upgrade magnon

conda-forge
-----------

If you prefer installation via the Conda package manager, you may use:

.. code-block:: bash

   conda install conda-forge::magnon

Source code (tarball)
---------------------

To install from a compressed source code file, start by downloading the latest tarball from [link]

Then extract the files, and create a symbolic link:

.. code-block::

   tar -xf [filename].tar.gz
   $ ln -s [filename] magnon

Thanks to the symbolic link, we can install using PyPI as:

.. code-block:: bash

   pip install /path/to/magnon

