import numpy as np
import numpy.linalg as npla
from numpy.typing import ArrayLike

import magnon
from magnon.common import default_distance_tol, normalize_scaled_coordinate



def cartesian_k_from_scaled(k_scaled, reciprocal_lattice_vectors):
    """
    Convert scaled k-vector(s) to Cartesian k-vector(s) using reciprocal lattice vectors.

    
    Parameters
    ----------
        k_scaled : ndarray
            Array of one or more k-vectors in scaled coordinates
        reciprocal_lattice_vectors : ndarray
            Reciprocal lattice vectors, as rows

        
    Returns
    ----------
        Array of k-vectors in Cartesian coordinates

    """

    k_scaled = np.atleast_2d(k_scaled).astype(np.float64)
    recip = np.asarray(reciprocal_lattice_vectors, dtype=np.float64)
    k_cartesian = k_scaled @ recip
    return k_cartesian[0] if k_cartesian.shape[0] == 1 and k_scaled.ndim == 1 else k_cartesian


def scaled_k_from_cartesian(k_cartesian, recip_lattice_vectors):
    """
    Given one or more Cartesian k vector(s), return the scaled k vector(s)

    
    Parameters
    ----------
        k_cartesian : ndarray
            Array of one or more k-vectors in Cartesian coordinates
        recip_lattice_vectors : ndarray
            Reciprocal lattice vectors, as rows

        
    Returns
    ----------
        Array of k-vectors in scaled coordinates

    """

    k_dir = np.atleast_2d(k_cartesian).astype(np.float64)
    recip = np.asarray(recip_lattice_vectors, dtype=np.float64)
    k_red = k_dir @ npla.inv(recip)
    return k_red[0] if k_red.shape[0] == 1 and k_dir.ndim == 1 else k_red

def reciprocal_lattice(lattice):
    """
    Compute reciprocal lattice basis vectors from lattice vectors

    
    Parameters
    ----------
        lattice : ndarray
            Array of lattice basis vectors, where each row is a basis vector

        
    Returns
    ----------
        recip_lattice : ndarray
            Reciprocal lattice basis vectors, where each row is a basis vector

    """

    def reciprocal_unit_vectors_3D(a1, a2, a3):
        cross = np.cross(a2, a3)
        triple = np.dot(a1, cross)
        return 2 * np.pi * cross / triple

    b1 = reciprocal_unit_vectors_3D(lattice[0], lattice[1], lattice[2])
    b2 = reciprocal_unit_vectors_3D(lattice[1], lattice[2], lattice[0])
    b3 = reciprocal_unit_vectors_3D(lattice[2], lattice[0], lattice[1])
    recip_lattice = [b1,b2,b3]

    recip_lattice = np.array(recip_lattice)

    return recip_lattice

def unit_spin_model_factors(spin_quantum_numbers,spins_are_unit):
    """
    Compute the Hamiltonian conversion factors for couplings from the unit spin model, or return 1s if model not in use

    
    Parameters
    ----------
        spin_quantum_numbers : ndarray
            1D array of spin quantum numbers
        spins_are_unit : bool
            Whether the unit spin model is being used

        
    Returns
    ----------
        spin_quantum_numbers : ndarray
            1D array of the spin quantum numbers to use under the chosen model
        factors : ndarray
            1D array of factors to apply to each element of the Hamiltonian for compatibility with the chosen model

        
    Notes
    ----------
        Some codes use unit spins rather than dimensionful spins whose magnitude is determined by the magnetic moment. The
        choice of model changes the definition of the exchange coupling constants. To enable compatibility with exchange
        coupling values from dimensionless spin models, this function may be used to apply an elementwise conversion along
        both axes of a linear spin wave Hamiltonian. The spin quantum numbers are set to 1, and the conversion along the
        atom index axis is returned.

    """
    
    if spins_are_unit:
        factors = 2 / np.where(spin_quantum_numbers != 0, spin_quantum_numbers, np.inf)
        spin_quantum_numbers = np.ones_like(spin_quantum_numbers)
    else:
        factors = np.ones_like(spin_quantum_numbers)

    return spin_quantum_numbers, factors
