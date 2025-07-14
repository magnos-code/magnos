import numpy as np
import spglib
from ase import Atoms
from scipy.spatial.distance import cdist

import magnon.interactions
from magnon.common import normalize_scaled_coordinate, reciprocal_lattice, default_distance_tol

def cell_transformation(original_cell, new_cell):
    """
    Computes the matrix transformation mapping from scaled coordinates in the original cell
    to scaled coordinates in the new cell.

    
    Parameters
    ----------
    original_cell : numpy.ndarray, dtype float
        The lattice vectors of the original cell with shape (3, 3).
    new_cell : numpy.ndarray, dtype float
        The lattice vectors of the new cell with shape (3, 3).

    
    Returns
    ----------
    transformation_matrix : numpy.ndarray, dtype float
        The matrix transforming a point from the old to the new basis with shape (3, 3).

    
    See Also
    ----------
    magnon.build.build_primitive_cell : Uses to map from conventional to primitive cell.

    
    Notes
    ----------

    .. math:: NTx = Ox

    for transformation_matrix :math:`T`, scaled coordinates :math:`x`, and old and new cell bases :math:`O,N` respectively.
    """
    return np.linalg.solve(new_cell.T, original_cell.T)

def build_primitive_cell(ase_atoms, interactions, distance_tol: float=default_distance_tol, no_idealize: bool = False):
    """
    From the provided Atoms object, attempts to find a primitive cell and rewrites the interactions for this primitive cell.

    Parameters
    ----------
        ase_atoms : ase.Atoms
            The atoms object from which to build the primitive cell.
        interactions : magnon.InteractionsList
            The InteractionList object for the original cell exchange couplings.
        distance_tol : float, optional
            The precision to use in determining equivalent positions.
        no_idealize : bool
            Whether Spglib should idealize the primitive cell.

    Raises
    ----------
        RuntimeError
            If the generated primitive cell includes atom properties not compatible with those of the retained atoms from the
            original cell.

    Returns
    ----------
        ase_prim_atoms : ase.Atoms object
            The ASE Atoms object for the primitive cell
        interactions : magnon.InteractionsList
            The InteractionList object for the primitive cell.
        transformation_matrix : numpy.ndarray, dtype float
            The matrix transformation mapping from scaled coordinates with respect to the original cell, with shape (3, 3).

    Notes
    ----------
        This 'idealizes' the lattice, meaning the primitive cell may be rotated in cartesian space. We do this because
        otherwise we can end up with strange primitive lattice vectors (e.g. negative or left-handed), which then causes
        great difficulty in plotting correct k-space paths. We next calculate the transformation matrix between the two
        cells in the next step which we will need for manipulating the interaction vectors.

        The spglib `standardize_cell` function does not use magnetic information. We cannot use
        spglib.get_magnetic_symmetry_dataset because this has no way to idealize the cell. This is an issue, for example,
        with antiferromagnets where often the magnetic cell must be larger than non-magnetic cell. We work around the
        problem here by assigning all atoms with the same atomic number and same magnetic moment a unique number. Spglib
        only sees numbers as type labels, not as atomic numbers, so this is fine. To avoid collisions with atomic
        numbers we add 1000x the magnetic moment index. This also means we can trivially remove the artificial type
        and restore the atomic numbers for the primitive cell data returned by spglib.
    """

    unique_moments = np.unique(ase_atoms.get_initial_magnetic_moments(), axis=0)
    magnetic_types = []
    for number, magmom in zip(ase_atoms.get_atomic_numbers(), ase_atoms.get_initial_magnetic_moments()):
        if np.allclose(magmom, 0):
            magnetic_types.append(number)
        else:
            magmom_index = np.where((unique_moments == magmom).all(axis=1))[0][0] + 1
            magnetic_types.append(number + 1000*magmom_index)

    spg_prim_cell = spglib.standardize_cell((ase_atoms.get_cell(), ase_atoms.get_scaled_positions(), magnetic_types), to_primitive=True, no_idealize=no_idealize)

    prim_lattice, prim_positions, prim_numbers = spg_prim_cell
    prim_numbers = prim_numbers % 1000

    # Calculate the mapping from the input cell to the output cell
    transformation_matrix = cell_transformation(ase_atoms.get_cell(), prim_lattice)

    # Calculate the mapping from the primitive cell to the original cell
    mapped_positions_scaled = np.dot(prim_positions, np.transpose(np.linalg.inv(transformation_matrix)))
    mapped_positions_scaled = normalize_scaled_coordinate(mapped_positions_scaled)

    distances = cdist(np.atleast_2d(mapped_positions_scaled), np.atleast_2d(ase_atoms.get_scaled_positions()))
    mapping = np.argmin(distances, axis=1)

    # Verify that the distances between the mapped and matched atoms are small
    matched_distance_indices = np.concatenate((np.atleast_2d(np.arange(distances.shape[0])), np.atleast_2d(mapping)), axis=0).transpose()
    matched_distances = distances[matched_distance_indices[:,0], matched_distance_indices[:,1]]
    if np.any(matched_distances > distance_tol):
        failing_mapped_atom = np.argwhere(matched_distances > distance_tol)
        raise RuntimeError(
            f"One or more matched atoms is not sufficiently close to its corresponding mapped atom:\n",
            f"Primitive cell atom was {failing_mapped_atom} at {np.atleast_2d(mapped_positions_scaled)[failing_mapped_atom].flatten()}\n",
            f"Original cell atom was {matched_distance_indices[failing_mapped_atom,1]} at "
            f"{np.atleast_2d(ase_atoms.get_scaled_positions())[matched_distance_indices[failing_mapped_atom,1]].flatten()}",)

    # Find the magnetic moments for the primitive atoms using the mapping
    prim_magmoms = np.zeros((len(prim_positions), 3))
    for i, original_index in enumerate(mapping):
        prim_magmoms[i] = ase_atoms.get_initial_magnetic_moments()[original_index]

    # Print the mapping information and check for consistency of atom type and magnetic moment
    for i, original_index in enumerate(mapping):
        original_atomic_number = ase_atoms.get_atomic_numbers()[original_index]
        original_magnetic_moment = ase_atoms.get_initial_magnetic_moments()[original_index]

        if prim_numbers[i] != original_atomic_number:
            raise RuntimeError("Mapping from primitive atomic numbers to original atoms is inconsistent.")

        if not np.allclose(prim_magmoms[i], original_magnetic_moment):
            raise RuntimeError("Mapping from primitive magnetic moments to original atoms is inconsistent.")

    # Calculate the reverse mapping from the original cell to the reverse cell
    reverse_mapped_positions_frac = np.dot(ase_atoms.get_scaled_positions(), np.transpose(transformation_matrix))
    reverse_mapped_positions_frac = normalize_scaled_coordinate(reverse_mapped_positions_frac)

    distances = cdist(reverse_mapped_positions_frac, prim_positions)
    reverse_mapping = np.argmin(distances, axis=1)

    # Print the reverse mapping information and check for consistency of atom type and magnetic moment
    for i, prim_index in enumerate(reverse_mapping):
        original_atomic_number = ase_atoms.get_atomic_numbers()[i]
        original_magnetic_moment = ase_atoms.get_initial_magnetic_moments()[i]

        if prim_numbers[prim_index] != original_atomic_number:
            raise RuntimeError("Mapping from primitive atomic numbers to original atoms is inconsistent.")

        if not np.allclose(prim_magmoms[prim_index], original_magnetic_moment):
            raise RuntimeError("Mapping from primitive magnetic moments to original atoms is inconsistent.")


    ase_prim_atoms = Atoms(
        symbols=prim_numbers,
        scaled_positions=prim_positions,
        cell=prim_lattice,
        magmoms=prim_magmoms
    )

    # Here I need the reverse mapping which maps the index in the original cell to the primitive cell
    prim_interactions = magnon.interactions.InteractionList([], atoms=ase_prim_atoms)

    for i, i_original_index in enumerate(mapping):
        for _, j_original_index, r_ij_scaled, J_ij in interactions.get_interactions(i_original_index):
            j = reverse_mapping[j_original_index]
            transformed_r_ij_scaled = np.dot(r_ij_scaled, np.transpose(transformation_matrix))

            original_r_ij_cartesian_length = np.linalg.norm(ase_atoms.get_cell().T @ r_ij_scaled)
            primitive_r_ij_cartesian_length = np.linalg.norm(ase_prim_atoms.get_cell().T @ transformed_r_ij_scaled)

            if not np.isclose(original_r_ij_cartesian_length, primitive_r_ij_cartesian_length, atol=default_distance_tol):
                raise RuntimeError("Mapping from primitive interactions to original interactions is inconsistent.")

            prim_interactions.insert_interaction(i, j, transformed_r_ij_scaled, J_ij)

    return (ase_prim_atoms, prim_interactions, transformation_matrix)

def convert_kpts(kpts, original_cell, new_cell):
    """
    Convert k-point(s) from scaled coordinates of one cell to another.
    
    Parameters
    ----------
        kpts : numpy.ndarray, dtype float
            The k-point(s) in scaled coordinates of the original_cell as an array of shape (3,) or (:math:`N_k`, 3)
        original_cell : numpy.ndarray, dtype float
            The lattice vectors of the original cell with shape (3,3)
        new_cell : numpy.ndarray, dtype float
            The lattice vectors of the new cell with shape (3,3)
        
    Returns
    ----------
        k_scaled_new_cell : numpy.ndarray, dtype float
            The k-point(s) in scaled coordinates of new_cell with shape (3,) or (:math:`N_k`, 3)
        
    See Also
    ----------
        magnon.build.cell_transformation : Returns the transformation matrix for scaled coordinates between two cells
    """

    kpts_array = np.atleast_2d(kpts).astype(np.float64)
    was_1d_input = kpts_array.shape[0] == 1 and np.asarray(kpts).ndim == 1

    original_recip_cell = reciprocal_lattice(original_cell)
    new_recip_cell = reciprocal_lattice(new_cell)
    transformation_matrix = cell_transformation(original_recip_cell, new_recip_cell)
    k_scaled_new_cell = np.dot(kpts_array, transformation_matrix)

    return k_scaled_new_cell[0] if was_1d_input else k_scaled_new_cell


