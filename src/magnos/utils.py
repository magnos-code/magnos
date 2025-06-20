import numpy as np
from ase import Atoms

def is_array_like(obj, ndim):
    """
    Check if the input object is an array-like structure with a specific number of dimensions.

    Parameters
    ----------
    obj : object
        Object to be checked.
    ndim : int
        Number of dimensions expected.

    Returns
    -------
    bool
        True if the object is array-like with the specified number of dimensions, False otherwise.

    Raises
    ------
    Exception
        If the object cannot be converted to a NumPy array.
    """
    try:
        arr = np.asarray(obj)
        return arr.ndim == ndim
    except Exception:
        return False

def all_same_shape(*arrays):
    """
    Check whether all input arrays have the same shape.

    Parameters
    ----------
    arrays : array-like
        Variable number of arrays to compare.

    Returns
    -------
    bool
        True if all arrays have the same shape, False otherwise.
    """
    shapes = [np.shape(a) for a in arrays]
    return len(set(shapes)) == 1

def all_same_length(*arrays):
    """
    Check whether all input arrays have the same length.

    Parameters
    ----------
    arrays : array-like
        Variable number of arrays to compare.

    Returns
    -------
    bool
        True if all arrays have the same length, False otherwise.
    """
    lengths = [len(a) for a in arrays]
    return len(set(lengths)) == 1

def ensure_vector_magnetic_moments(atoms):
    """
    Ensure that the ASE Atoms object has vector (3D) magnetic moments for each atom.

    Parameters
    ----------
    atoms : ase.Atoms
        ASE atoms object to be validated.

    Returns
    -------
    atoms : ase.Atoms
        ASE Atoms object with vector magnetic moments.

    Notes
    -----
    ASE allows `initial_magnetic_moments` to be either scalar or vector. This function
    guarantees vector format, which is useful for handling non-collinear magnetism.
    """
    if atoms.get_initial_magnetic_moments().ndim == 1:
        vector_magnetic_moments = np.atleast_2d(
            np.einsum("i,j->ij", atoms.get_initial_magnetic_moments(), np.array([0, 0, 1]))
        )
        atoms = Atoms(
            atoms.get_chemical_symbols(),
            positions=atoms.get_positions(),
            cell=atoms.get_cell()
        )
        atoms.set_initial_magnetic_moments(vector_magnetic_moments)
    return atoms
