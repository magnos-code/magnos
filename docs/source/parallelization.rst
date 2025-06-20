Parallelization
===============

Calculating a bandstructure over a path of many points in reciprocal space is an embarrassingly parallel task, as the process of
calculating the bandstructure at one k-point does not depend on the result at another k-point.

The bandstructure calculation is parallelized using multiprocessing. This can be activated by setting :code:`num_threads` to
an integer other than :math:`1` when initializing the instance of :code:`MagnonSpectrum`.