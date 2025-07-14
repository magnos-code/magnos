Development
===========

Magnon is designed to work as a component of larger workflows. This underpins its development principles:

* Integrate with existing physics libraries for handling atomic structures (e.g. ASE)
* Mirror design principles of existing physics libraries (e.g. the MagnonSpectrum class has the same structure as ASE's Phonons class)
* Keep things clean and simple so that others can build on top of Magnon. This also mandates rigorous documentation, comprehensive tutorials and readable code

We hope you'll choose to extend our capabilities!

GitHub
----------

The source code for Magnon is hosted on `GitHub <https://github.com/>`_. This is the place to:

* Raise bug reports
* Suggest improvements/new features
* Request to merge your own extensions

Cloning the code
++++++++++++++++

To install directly from the source code, you may clone from our GitHub repository:

.. code-block:: bash

   git clone -b stable https://github.com/

Then either:

1. Install the package by ensuring the path to the :code:`magnon` directory is included in your system's :code:`PYTHONPATH`
variable:

.. code-block:: bash

   export PYTHONPATH="$PYTHONPATH:/path/to/magnon"

Or...

2. Install using PyPI from the downloaded source code:

.. code-block:: bash

   pip install /path/to/magnon

.. note::

   You can use :code:`pip` with :code:`-e` to install in editable mode, so that your changes are immediately available in the
   installed module.
