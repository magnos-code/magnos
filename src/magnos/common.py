import numpy as np
from ase.calculators.calculator import Parameters
from numpy.typing import ArrayLike

# Shared default tolerance (used across lattice, build, and interactions)
default_distance_tol = 1e-5
"""
Default precision used in determining equivalent distances, value = 1e-5
"""


def normalize_scaled_coordinate(r, distance_tolerance: float = default_distance_tol):
    """
    Normalizes a scaled coordinate vector or an array of scaled coordinates
    to the range [0, 1) for each component.

    Parameters
    ----------
        r : array_like
            A scaled coordinate vector or an array of scaled coordinates
        distance_tolerance : float
            Precision used for rounding near cell boundaries.

    Returns
    ----------
        numpy.ndarray, dtype = np.float64
            A normalized scaled coordinate vector or an array of normalized scaled coordinates
    """
    r = np.asarray(r, dtype=np.float64)
    orig_shape = r.shape
    r = np.atleast_2d(r)
    r = r - np.floor(r)
    r[np.isclose(r, 1.0, atol=distance_tolerance)] = 0.0
    return r[0] if orig_shape == (3,) else r


def modulo_lattice(x: ArrayLike, distance_tol: float = default_distance_tol) -> np.ndarray:
    """
    Reduces a lattice vector x modulo 1, ensuring that values numerically close to 1 or 0 are wrapped to 0.

    Parameters
    ----------
        x : array_like
            A lattice vector
        distance_tol : float
            Precision used for rounding near cell boundaries.

    Returns
    ----------
        numpy.ndarray
            A reduced lattice vector x modulo 1

    """
    arr = np.asarray(x) % 1.0
    arr[np.isclose(arr, 1.0, atol=distance_tol) | np.isclose(arr, 0.0, atol=distance_tol)] = 0.0
    return arr


def squared_distance(a, b):
    """
    Computes squared Euclidean distance between two points.

    Parameters
    ----------
        a : array_like
            Coordinates of the first point
        b : array_like
            Coordinates of the second point

    Returns
    -------
        float
            Squared Euclidean norm of distance vector :math:`a-b`.
    """

    c = a - b
    return np.dot(c, c)


def find_site_index(position, scaled_atom_positions, distance_tol: float = default_distance_tol):
    """
    Given a position in scaled coordinates, find the equivalent site index in the unit cell.

    Parameters
    ----------
        position : array_like
            A scaled coordinate vector to search for
        scaled_atom_positions : array_like
            An array of scaled coordinates in which to search
        distance_tol : float
            Precision used for determining equivalent positions

    Raises
    ----------
        ValueError
            If a matching site could not be found

    Returns
    ----------
        int
            The index of site which matches the coordinates
    """
    r = normalize_scaled_coordinate(position)
    for i, pos in enumerate(scaled_atom_positions):
        if np.allclose(r, pos, atol=distance_tol):
            return i
    raise ValueError("Could not find a matching site")


def lattice_translation_vector(r_ij_scaled: np.ndarray, r_i_scaled: np.ndarray, r_j_scaled: np.ndarray,
                               distance_tol: float = default_distance_tol) -> np.ndarray:
    """
    Returns the integer lattice translation vector T of an arbitrary vector r_scaled.

    Parameters
    ----------
        r_ij_scaled : array_like
            A scaled distance vector from position i to position j.
        r_i_scaled : array_like
            The scaled coordinates of point i
        r_j_scaled : array_like
            The scaled coordinates of point j
        distance_tol : float
            Precision used for checking if a value is an integer

    Raises
    ----------
        ValueError
            If the three vectors are inconsistent and do not result in an integer scaled lattice translation vector

    Returns
    ----------
        numpy.ndarray, dtype = int
            The integer lattice translation vector
    """
    q_ij = r_ij_scaled - (r_j_scaled - r_i_scaled)
    if not np.allclose(np.rint(q_ij), q_ij, atol=distance_tol):
        raise ValueError(f"lattice translation vector {q_ij} is non-integer")
    nearest_integers = np.round(q_ij)
    floored_values = np.floor(q_ij)
    mask = np.isclose(q_ij, nearest_integers, atol=distance_tol)
    translation_vector = np.where(mask, nearest_integers, floored_values).astype(int)
    return translation_vector


def reciprocal_lattice(lattice):
    """
    Compute reciprocal lattice basis vectors from lattice vectors.

    Parameters
    ----------
        lattice : array_like
            A 3x3 array of lattice basis vectors; each row is a vector

    Returns
    ----------
        numpy.ndarray
            A 3x3 array of reciprocal lattice basis vectors; each row is a vector
    """
    def reciprocal_unit_vectors_3D(a1, a2, a3):
        cross = np.cross(a2, a3)
        triple = np.dot(a1, cross)
        return 2 * np.pi * cross / triple

    b1 = reciprocal_unit_vectors_3D(lattice[0], lattice[1], lattice[2])
    b2 = reciprocal_unit_vectors_3D(lattice[1], lattice[2], lattice[0])
    b3 = reciprocal_unit_vectors_3D(lattice[2], lattice[0], lattice[1])
    return np.array([b1, b2, b3])

