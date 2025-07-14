Theory
===============

Magnon is a Python package for the quantitative simulation of magnons, i.e., quantized spin-wave excitations representing the coherent precession of spins around their equilibrium direction.

In crystalline magnetic solids, magnons can be described using the second-quantization formalism. In particular, magnon creation and annihilation operators are labeled by: (i) a wavevector that characterizes the periodicity of the spin wave; and (ii) a band index, which distinguishes between different excitation patterns in crystals with multiple magnetic sites in the primitive cell.

Suggested reading
-----------------
To follow the theory in this section, readers should have a background in quantum mechanics, second quantization, and condensed matter physics. Some useful references for these topics are listed below:

1. Quantum mechanical operators: B.H. Bransden and C.J. Joachain, "Quantum Mechanics", 2ed, Pearson (2000), chapter 5
2. Second quantization: J.M. Ziman, "Electrons and Phonons", Oxford University Press (2001), chapter 1
3. Crystal and reciprocal lattices: C. Kittel, "Introduction to Solid State Physics", 8ed, Wiley (2005), chapters 1 and 2

The magnon bandstructure can be derived using Linear Spin Wave Theory (LSWT) at the quantum level, or semiclassical spin dynamics - Magnon uses the former.
For an introduction to magnons, see the derivation on pages 330-333 of Kittel, referenced above.

Organization of the theory
---------------------------

In the 'Background' section, we summarize the origins of angular momentum in atoms, the connection to magnetic moments, the quantum mechanical operators for spin, and the Heisenberg Hamiltonian for spins.

In the 'Linear Spin Wave Theory' section, we derive the general form of the Heisenberg Hamiltonian used by Magnon, and discuss all the steps necessary to diagonalize it and obtain the magnon bandstructure.

After discussing the basic theoretical details, we will present several hands-on tutorials.

.. toctree::
   :maxdepth: 1

   background
   linear-spin-wave-theory
