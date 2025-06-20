Theory
===============

Magnos is a Python package for the quantitative study of magnons. Spin waves are elementary excitations of
a spin system, representing oscillations in the relative orientations of magnetic moments. Magnons are the quasiparticle
representation of these spin waves. This means that the excitation of a spin wave on top of some ground-state
spin configuration can be alternatively represented as the creation of a magnon.

The creation of a magnon costs some energy, the amount depending on the wavevector and the particular band. The wavevector
describes the periodicity of the spin wave compared to the lattice, and the band determines the particular eigenstate at that
wavevector.

Suggested reading
-----------------
To follow the theory in this section, the necessary background is courses in quantum mechanics, second quantization and
condensed matter physics. Some useful references for these areas are:

1. Quantum mechanical operators: B.H. Bransden and C.J. Joachain, "Quantum Mechanics", 2ed, Pearson (2000), chapter 5
2. Second quantization: J.M. Ziman, "Electrons and Phonons", Oxford University Press (2001), chapter 1
3. Crystal and reciprocal lattices: C. Kittel, "Introduction to Solid State Physics", 8ed, Wiley (2005), chapters 1 and 2

The magnon bandstructure can be derived using Linear Spin Wave theory at the quantum level, or semiclassical spin dynamics - Magnos uses the former.
For an introduction to magnons, see the derivation using the latter on pages 330-333 of Kittel, referenced above.

In this section
---------------

In the 'Background' section, we recap the origins of angular momentum in atoms, the connection to magnetic moments, the quantum mechanical
operators for spin, and the theoretical justification for the Heisenberg Hamiltonian for spins.

In the 'Linear Spin Wave Theory' section, we derive the Hamiltonian used by Magnos to obtain the magnon bandstructure.

After familiarizing yourself with the details in this section, you'll be ready to try our hands-on tutorials.

.. toctree::
   :maxdepth: 1

   background
   linear-spin-wave-theory
