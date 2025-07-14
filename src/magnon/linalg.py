import numpy as np
from numpy.typing import ArrayLike

import magnon


def normalised_vector(vec: ArrayLike) -> np.ndarray:
    """
    Return a normalised vector in the direction of `vec`.
    
    Parameters
    ----------
        vec : array_like
            The vector to normalize.
    
    Returns
    ----------
        out : ndarray
            Normalised vector with the same direction as `vec`. If the input norm is zero, returns a vector of zeros with the same shape as `vec`.

    Examples
    --------
    >>> normalised_vector([3, 4])
    array([0.6, 0.8])
    >>> normalised_vector([0, 0])
    array([0., 0.])
    """
    arr = np.asarray(vec)
    norm = np.linalg.norm(arr)
    if norm == 0:
        return np.zeros_like(arr)
    return arr / norm


def generate_orthonormal_vectors(vec: ArrayLike) -> tuple[np.ndarray, ...]:
    """
    Generate an orthonormal set of vectors including the given input vector.
    
    Given a 2D or 3D input vector, this function returns a set of orthonormal
    vectors where the first vector is the normalised input vector, and the
    remaining vectors are mutually orthogonal and normalised.
    
    
    Parameters
    ----------
        vec : array_like
            Input vector of length 2 or 3.
    
    
    Returns
    ----------
        orthonormals : tuple of ndarray
            Tuple containing orthonormal vectors:

            For 2D input: (v0, v1), where v0 is the normalised input vector,
            and v1 is perpendicular to v0.

            For 3D input: (v0, v1, v2), where v0 is the normalised input vector,
            and v1, v2 are mutually orthogonal and normalised.
    
    
    Raises
    ----------
        AssertionError
            If the length of the input vector is not 2 or 3.
    
    Examples
    --------
    >>> generate_orthonormal_vectors([1, 0, 0])
    (array([1., 0., 0.]), array([0., 0., 1.]), array([0., 1., 0.]))
    
    >>> generate_orthonormal_vectors([0, 1])
    (array([0., 1.]), array([1., 0.]))
    """
    assert len(vec) == 2 or len(vec) == 3

    rng = np.random.default_rng()

    # make a normalised copy so that we don't change the original vec
    normalised_vec = normalised_vector(vec)

    # guess a non-collinear vector, keep trying until one is found
    non_collinear_vec = rng.standard_normal(normalised_vec.shape)

    while np.all(np.cross(non_collinear_vec, normalised_vec) == 0.0):
        non_collinear_vec = rng.standard_normal(normalised_vec.shape)

    # use cross product to get perpendicular vectors
    perp_vec_1 = normalised_vector(np.cross(non_collinear_vec, normalised_vec))

    if len(vec) == 3:
        perp_vec_2 = normalised_vector(np.cross(normalised_vec, perp_vec_1))
        return normalised_vec, perp_vec_1, perp_vec_2
    else:
        return normalised_vec, perp_vec_1

def rotation_matrix_pair(axis: ArrayLike, angle: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Return the generalized rotation matrix describing a rotation of [angle] radians around [axis], and a rotation of -[angle] radians around [axis]

    This is an implementation of Rodrigues' rotation formula [1]_.

    Parameters
    --------
        axis : array_like
            The axis around which the rotation is to occur
        angle: float
            The angle of the rotation in radians
    
    Returns
    --------
        R : numpy.ndarray
            The rotation matrix
    
    References
    ----------
    
    .. [1] https://arxiv.org/pdf/1810.02999
    
    """
    axis = np.asarray(axis, dtype=np.float64)
    axis = axis / np.linalg.norm(axis)
    c = np.cos(angle)
    s = np.sin(angle)
    t = 1 - c
    x, y, z = axis

    # Rodrigues' rotation formula
    K = np.array([
        [0, -z, y],
        [z, 0, -x],
        [-y, x, 0]
    ])
    outer = np.outer(axis, axis)

    def rot(sign=1):
        return c * np.eye(3) + sign * s * K + t * outer

    R_forward = rot(sign=1)
    R_reverse = rot(sign=-1)
    return R_forward, R_reverse

def paraidentity(size: int, dtype: np.dtype=np.float64) -> np.ndarray:
    """
    Create a paraidentity (:math:`\\sigma_3`) matrix of shape (2*size, 2*size).
    
    The paraidentity matrix, also known as the generalized Pauli :math:`\\sigma_3` matrix,
    is a diagonal matrix with the first `size` diagonal entries set to 1, and
    the next `size` diagonal entries set to -1.
    
    
    Parameters
    ----------
        size : int
            The block size. The output matrix will have shape (2*size, 2*size).
        dtype : data-type, optional
            Desired output data type (default is `np.float64`).
    
    
    Returns
    ----------
        out : ndarray
            Paraidentity matrix of shape (2*size, 2*size) and specified dtype.
    
    Examples
    --------
    >>> paraidentity(2)
    array([[ 1.,  0.,  0.,  0.],
    [ 0.,  1.,  0.,  0.],
    [ 0.,  0., -1.,  0.],
    [ 0.,  0.,  0., -1.]])
    """

    return np.diag(np.array([1]*size + [-1]*size, dtype=dtype))


def is_purely_real(arr: ArrayLike, zero_tol: float=magnon.default_numerical_tol):
    """
    Test whether all elements of an array are real within a tolerance.
    
    
    Parameters
    ----------
        arr : array_like
            Input array to be tested.
        zero_tol : float, optional
            Tolerance for considering the imaginary part as zero (default is 1e-8).
    
    
    Returns
    ----------
        is_real : bool
            True if all elements have imaginary part less than `tol` in magnitude, False otherwise.
    
    Examples
    --------
    >>> import numpy as np
    >>> is_purely_real([1, 2, 3])
    True
    >>> is_purely_real([1+1e-9j, 2, 3])
    True
    >>> is_purely_real([1+1e-6j, 2, 3], zero_tol=1e-7)
    False
    >>> is_purely_real([1+0j, 2+0j, 3+0j])
    True
    """
    return not np.any(np.abs(np.imag(np.asarray(arr))) > zero_tol)