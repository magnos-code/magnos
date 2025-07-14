import numpy as np


def band_structures_are_equal(kpts1, energies1, kpts2, energies2, tol=1e-8):
    """
    Compares two sets of k-points and energies to check if they are equivalent
    (within a tolerance for each point) allowing for reordering of k-points.
    Energies along the M dimension are assumed to be sorted.

    Parameters:
        kpts1: np.ndarray
            First set of k-points (N, 3 array).
        energies1: np.ndarray
            First set of energies (1, N, M array).
        kpts2: np.ndarray
            Second set of k-points (N, 3 array).
        energies2: np.ndarray
            Second set of energies (1, N, M array).
        tol: float
            Tolerance for numerical differences.

    Returns:
        bool: True if the two sets of points (kpts and energies) are equivalent within tolerance, False otherwise.
    """
    # Ensure inputs are NumPy arrays
    kpts1 = np.asarray(kpts1)
    energies1 = np.asarray(energies1).squeeze(axis=0)  # Shape becomes (N, M)
    kpts2 = np.asarray(kpts2)
    energies2 = np.asarray(energies2).squeeze(axis=0)  # Shape becomes (N, M)

    # Check that shapes are compatible
    if kpts1.shape != kpts2.shape or energies1.shape != energies2.shape:
        return False

    # Combine k-points and energies for sorting
    points1 = np.hstack((kpts1, energies1))  # Shape becomes (N, 3 + M)
    points2 = np.hstack((kpts2, energies2))  # Shape becomes (N, 3 + M)

    # Sort both sets of points by k-points (assumes sorting along energies is not required)
    sorted1 = points1[np.lexsort(points1.T[:3][::-1])]  # Sort by the first three columns (k-points)
    sorted2 = points2[np.lexsort(points2.T[:3][::-1])]  # Sort by the first three columns (k-points)

    # Compare using numpy.isclose
    return np.all(np.isclose(sorted1, sorted2, atol=tol))
