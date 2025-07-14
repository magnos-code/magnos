Tips
=============

Reproducing literature bandstructures
-------------------------------------

It can be difficult to reproduce magnon bandstructures in research papers. These steps can help overcome some common pitfalls:

1. Check the Hamiltonian convention - if the difference is a factor of 0.5 or 2, it may be that the reference uses a different prefactor for the Hamiltonian. By default, Magnon sums over all pairs of site indices :math:`i` and :math:`j`, which corresponds to the Heisenberg Hamiltonian multiplied by :math:`2`. You can change the prefactor by setting :code:`ham_prefactor` when initializing an instance of :code:`MagnonSpectrum`.
2. Check whether the exchange coupling values are computed using a model with unit vector (dimensionless) spins. By default, Magnon assumes that they are. If not, you can set :code:`spins_are_unit = False` when initializing the :code:`MagnonSpectrum`.
3. Check that the exchange coupling values are not in an unexpected unit, such as mRy rather than meV.
4. Check that the minimal set of couplings you supplied were sufficient to generate all of the couplings expected due to symmetry. Specify more couplings if needs be.
5. Check that you are using the correct path in reciprocal space

   a. Is the path used in the reference for a conventional or primitive cell?
   b. Do they use different symbols to those used elsewhere?
   c. Do they use Cartesian or scaled points?

Finding structure files
-----------------------

See the `International Crystal Structure Database <https://icsd.products.fiz-karlsruhe.de/>`_ for CIF files for many structures referenced in research literature.

