Magnon: Magnon Bandstructures in Python
=======================================

**Magnon** is a Python library for calculating magnon bandstructures from exchange coupling interactions. It is designed
to be intuitive and easy-to-use so that you can spend less time *calculating* bandstructures and more time *using* them
to do interesting science.

Install
-------

.. code-block::

   pip install magnon

Bandstructures in under 10 lines of code
----------------------------------------

It should be extremely easy to take exchange coupling data from the literature and reproduce the magnon bandstructure. Using
Magnon, the simple script below can be used to obtain the bandstructure of BCC iron in the primitive cell..

.. code-block::

   import magnon
   import magnon.build

   atoms, interactions = magnon.input.create_interacting_system('Iron_BCC.poscar', 'Iron_BCC_mag_moments', 'Iron_BCC_exchange', 2)

   interactions = interactions.symmetrize(atoms)
   primitive_atoms, primitive_interactions, _ = magnon.build.build_primitive_cell(atoms, interactions)

   magnon = magnon.MagnonSpectrum(primitive_atoms,primitive_interactions,num_threads=4)

   path = primitive_atoms.get_cell().bandpath(path='GNPGHN', npoints=180)

   bstruct = magnon.get_band_structure(path)
   bstruct.plot(emin=0, emax=0.7, filename='Iron_BCC_bands.png')

.. figure:: Iron_BCC_bandstructure.png

   The bandstructure of BCC Iron calculated using the above example script.

Why use Magnon?
---------------

Our core development principles are:

* **Easy to use** - Magnon bandstructures can be computed using less than 10 lines of Python code!
* **Scalable** - For handling complex cells, symmetries, and couplings.
* **Integrated** - Interfaces with tools such as ASE to accelerate scientific workflows. The MagnonSpectrum class mirrors the structure of ASE's Phonons class.
* **Adaptable** - Allows conventions for working with dimensionless or dimensionful magnetic moments, different Hamiltonian conventions, many input formats.
* **Open** - Built for collaboration.

The code offers:

* ASE-enabled reading of many widely used structure file formats
* Ability to read complex exchange interaction data, and rigorously check symmetry requirements
* Option to easily augment non-magnetic structure files with a magnetic moment data file
* Rapid set-up routine by simultaneously reading structural, magnetic and coupling data
* Powerful InteractionList class for working efficiently with exchange coupling data at scale
* In-built symmetrisation routines for generating full exchange couplings from a minimal description
* Unit cell standardisation to convert couplings to those of the primitive cell
* Supports ferromagnetic, antiferromagnetic and non-collinear spin configurations
* Adapts to different conventions for the Heisenberg Hamiltonian prefactor
* Able to convert from exchange couplings from unit vector spin models
* Integrated MagnonSpectrum class for easy generation of reciprocal space paths and their bandstructures. Built with the same structure as the ASE Phonons class.

