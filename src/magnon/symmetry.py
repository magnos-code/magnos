from itertools import product
from typing import Any
from collections.abc import Collection

import ase
import numpy as np
import spglib
from numpy.typing import ArrayLike

import magnon
from magnon.common import modulo_lattice, squared_distance, normalize_scaled_coordinate

class SymOp:
    """
    Represents a symmetry operation in three-dimensional space, consisting of a rotation and a translation.

    Parameters
    ----------
        rotation : array_like, shape (3, 3)
            Rotation matrix.
        translation : array_like, shape (3,)
            Translation vector.

    Notes
    ----------

    A symmetry operation is defined as a pair :math:`(R, \\mathbf{t})` where:

    * :math:`R` is a 3x3 rotation matrix
    * :math:`\\mathbf{t}` is a 3-element translation vector

    This class supports standard operations such as composition (multiplication), inversion,
    application to coordinates, and introspection of geometric properties such as rotation axes and fixed points.
    """
    def __init__(self, rotation: ArrayLike, translation: ArrayLike, symmetry_tol: float=magnon.default_numerical_tol):
        """
        Initialize a symmetry operation.


        """
        rotation = np.asarray(rotation, dtype=float)
        translation = np.asarray(translation, dtype=float)
        assert rotation.shape == (3, 3), "rotation must be 3x3"
        assert translation.shape == (3,), "translation must be length 3"
        self.rotation = rotation
        self.translation = translation
        self.symmetry_tol = symmetry_tol

    def __mul__(self, other: 'SymOp') -> 'SymOp':
        """
        Compose two symmetry operations via multiplication.

        Parameters
        ----------
        other : SymOp
            The symmetry operation to compose on the right.

        Returns
        ----------
        SymOp
            The composed symmetry operation.

        Notes
        ----------
        For two symmetry operations,

        .. math::

           S_1 x = R_1 x + T_1

        .. math::

           S_2 x = R_2 x + T_2

        the action of the product

        .. math::

           S_1 S_2 x = S_1 (R_2 x + T_2)

        .. math::

           = R_1 (R_2 x + T_2) + T_1

        .. math::

           = R_1 R_2 x + R_1 T_2 + T_1

        so that the composite symmetry operations consists of rotation :math:`R_1 R_2` and translation :math:`R_1 T_2 + T_1`.

        """
        return SymOp(self.rotation @ other.rotation,
                     self.translation + self.rotation @ other.translation)

    def __matmul__(self, pos: ArrayLike) -> np.ndarray:
        """
        Apply this symmetry operation to a coordinate vector.

        Parameters
        ----------
        pos : array_like, shape (3,)
            The position vector to transform.

        Returns
        -------
        np.ndarray
            The transformed coordinate.
        """
        return self.rotation @ np.asarray(pos) + self.translation

    def __add__(self, translation: ArrayLike) -> 'SymOp':
        """
        Return a copy of the symmetry operation with the translation vector incremented.

        Parameters
        ----------
        translation : array_like, shape (3,)
            Translation vector to add.

        Returns
        -------
        SymOp
            The symmetry operation with updated translation.
        """
        return SymOp(self.rotation, self.translation + np.asarray(translation))

    def __sub__(self, translation: ArrayLike) -> 'SymOp':
        """
        Return a copy of the symmetry operation with the translation vector decremented.

        Parameters
        ----------
        translation : array_like, shape (3,)
            Translation vector to subtract.

        Returns
        -------
        SymOp
            The symmetry operation with updated translation.
        """
        return SymOp(self.rotation, self.translation - np.asarray(translation))

    def __pow__(self, exponent: int) -> 'SymOp':
        """
        Return the symmetry operation raised to the given non-negative integer power.

        Parameters
        ----------
        exponent : int
            The exponent (must be non-negative).

        Returns
        -------
        SymOp
            The composed operation applied `exponent` times.
        """
        if exponent < 0:
            raise ValueError("exponent must be ≥ 0")
        if exponent == 0:
            return SymOp.identity()
        result = self
        for _ in range(exponent - 1):
            result = result * self  # Use non-mutating multiply
        return result

    def __eq__(self, other: Any) -> bool:
        """
        Check equality of two symmetry operations with numerical tolerance.

        Parameters
        ----------
        other : Any
            The object to compare with.

        Returns
        -------
        bool
            True if the operations are equal within tolerance.
        """
        return (isinstance(other, SymOp) and
                np.allclose(self.rotation, other.rotation, atol=self.symmetry_tol) and
                np.allclose(self.translation, other.translation, atol=self.symmetry_tol))

    def __hash__(self) -> hash:
        """
        Compute a hash value based on the rounded rotation and translation components.

        Returns
        -------
        int
            Hash value.
        """
        return hash((
            tuple(np.round(self.rotation.flatten(), 8)),
            tuple(np.round(self.translation, 8))
        ))

    def __repr__(self) -> str:
        """
        Return an unambiguous string representation.

        Returns
        -------
        str
            String representation of the symmetry operation.
        """
        return f"SymOp(rotation={self.rotation.tolist()},translation={self.translation.tolist()})"

    def __str__(self) -> str:
        """
        Return a formatted string showing the rotation matrix and translation vector.

        Returns
        -------
        str
            Human-readable representation of the symmetry operation.
        """
        return "\n".join(
            " ".join(f"{x: 3.3f}" for x in row) + f"  | {t: 3.3f}"
            for row, t in zip(self.rotation, self.translation)
        )

    @classmethod
    def identity(cls) -> 'SymOp':
        """
        Return the identity symmetry operation.

        Returns
        -------
        SymOp
            The identity operation (rotation = I, translation = 0).
        """
        return cls(np.eye(3), np.zeros(3))

    def inverse(self) -> 'SymOp':
        """
        Compute the inverse of this symmetry operation.

        Returns
        -------
        SymOp
            The inverse operation :math:`(R^{-1}, -R^{-1}t)`.
        """
        rot_inv = self.rotation.T
        return SymOp(rot_inv, -rot_inv @ self.translation)

    def rotation_order(self) -> int:
        """
        Return the order of the rotational part of this symmetry operation.

        The order is the smallest positive integer n such that :math:`R^n = I`.

        Returns
        ----------
        int
            The order of the rotation.

        Notes
        ----------
        If a rotation is combined with inversion (i.e. an improper rotation), it will have a determinant of -1.
        In this case, the order of the rotational part (excluding the inversion) is returned as a negative integer.
        If this order is odd, it must be doubled to ensure that the inversion part when raised to the exponent gives
        the identity.
        """
        det = round(np.linalg.det(self.rotation))
        tr = round(np.trace(self.rotation))
        if det == -1:
            return {-3: -1, -2: -6, -1: -4, 0: -3, 1: -2}.get(tr, None)
        elif det == 1:
            return {-1: 2, 0: 3, 1: 4, 3: 1}.get(tr, None)
        raise ValueError('Unknown symmetry operation: trace={}, det={}'.format(tr, det))

    def intrinsic_translation(self) -> np.ndarray:
        """
        Compute the intrinsic translation vector of this symmetry operation.

        This represents the translational component that remains after applying the operation repeatedly.

        Returns
        -------
        np.ndarray
            The intrinsic translation vector.
        """
        order = self.rotation_order()
        n = -2 * order if order in (-1, -3) else abs(order)
        symop_n = self ** n
        return symop_n.translation / n

    def fixed_point(self) -> np.ndarray:
        """
        Return a point that is invariant under this symmetry operation.

        Returns
        -------
        np.ndarray
            A fixed point in space under the operation.
        """
        return self.translation - self.intrinsic_translation()

    def rotation_axis(self, tol: float=magnon.default_numerical_tol) -> np.ndarray | None:
        """
        Return the axis of rotation as a unit vector, if defined.

        Parameters
        ----------
        tol : float
            Tolerance for identifying unit eigenvalues.

        Returns
        -------
        np.ndarray or None
            A unit vector along the rotation axis, or None if not well-defined.
        """
        order = self.rotation_order()
        if not abs(order) > 1:
            return None

        eigvals, eigvecs = np.linalg.eig(self.rotation)
        # Look for eigenvalue close to 1
        idx = np.where(np.isclose(eigvals, 1, atol=tol))[0]
        if len(idx) == 0:
            return None  # No unique axis
        axis = np.real(eigvecs[:, idx[0]])
        axis = axis / np.linalg.norm(axis)
        # Canonical orientation: z, then y, then x > 0 (Grosse-Kunstleve convention)
        for idx in (2, 1, 0):  # z, y, x
            val = axis[idx]
            if abs(val) > tol:
                if val < 0:
                    axis = -axis
                break
        return axis

    def is_translation_only(self, tol: float=magnon.default_numerical_tol) -> bool:
        """
        Check whether the symmetry operation is a pure translation.

        Parameters
        ----------
        tol : float
            Tolerance for determining if the rotational part is close to the identity.

        Returns
        -------
        bool
            True if the rotation is the identity and the translation is non-zero modulo lattice.
        """
        if not np.isclose(np.trace(self.rotation), 3, atol=tol):
            return False
        if not np.allclose(self.rotation, np.eye(3), atol=tol):
            return False

        translation_mod = normalize_scaled_coordinate(self.translation)
        return np.dot(translation_mod, translation_mod) > tol ** 2

    def is_identity(self, tol: float=magnon.default_numerical_tol) -> bool:
        """
        Check whether the symmetry operation is the identity.

        Parameters
        ----------
        tol : float
            Tolerance for determining if the operation is close to the identity.

        Returns
        -------
        bool
            True if the operation is (I, 0).
        """
        # Fast cheap tests first:
        if np.any(np.abs(self.translation) > tol):
            return False
        if not np.isclose(np.trace(self.rotation), 3, atol=tol):
            return False
        if not np.allclose(self.rotation, np.eye(3), atol=tol):
            return False
        return True

    def as_homogeneous_matrix(self) -> np.ndarray:
        """
        Return the 4x4 homogeneous transformation matrix for this symmetry operation.

        Returns
        -------
        np.ndarray
            A 4x4 matrix combining rotation and translation.
        """
        return np.block([
            [self.rotation, self.translation.reshape(3, 1)],
            [np.zeros((1, 3)), np.ones((1, 1))]
        ])


def symop_product_modulo_lattice(a: SymOp, b: SymOp) -> SymOp:
    """
    Multiply two symmetry operations with translation wrapped modulo lattice.

    This function performs composition of two symmetry operations of the form
    :math:`(R_a, \\mathbf{t}_a)` and :math:`(R_b, \\mathbf{t}_b)`, and reduces the
    resulting translation vector modulo 1. The operation is defined as:

    .. math::

       R = R_a R_b

    .. math::

       \\mathbf{t} = (\\mathbf{t}_a + R_a \\mathbf{t}_b) \\; \\text{mod} \\; 1

    This is useful when working with symmetry operations in scaled coordinates
    where translations differing by a lattice vector are equivalent.

    Parameters
    ----------
    a : SymOp
        The first symmetry operation.

    b : SymOp
        The second symmetry operation.

    Returns
    -------
    SymOp
        The resulting symmetry operation after composition, with translation modulo 1.

    See Also
    ----------
        `magnon.symmetry.SymOp.__mul__`
    """
    rot = a.rotation @ b.rotation
    trans = normalize_scaled_coordinate(a.translation + a.rotation @ b.translation)
    return SymOp(rot, trans)


def spgcell(ase_atoms, magnetic=False) -> (np.ndarray, np.ndarray, np.ndarray, np.ndarray | None):
    """
    Convert an ASE Atoms object into a spglib-compatible tuple representation.

    This function returns the standard input format required by spglib:
    (lattice, positions, atomic numbers [, magnetic moments]), depending on
    whether magnetic symmetry detection is requested.

    Parameters
    ----------
    ase_atoms : ase.Atoms
        An ASE Atoms object representing the crystal structure.

    magnetic : bool, optional
        If True, include initial magnetic moments in the output tuple for magnetic symmetry analysis.
        Default is False.

    Returns
    -------
    tuple
        A tuple containing:

        * lattice : np.ndarray, shape (3, 3) - Lattice vectors as rows of a matrix.

        * positions : np.ndarray, shape (N, 3) - Scaled atomic positions.

        * numbers : np.ndarray, shape (N,) - Atomic numbers.
        
        * moments : np.ndarray, shape (N,), optional - Magnetic moments, only returned if `magnetic=True`.

    Notes
    -----
    This function is useful for preparing input for `spglib.get_symmetry`, `get_spacegroup`, etc.
    """

    if magnetic:
        return (
            ase_atoms.get_cell(),
            ase_atoms.get_scaled_positions(),
            ase_atoms.get_atomic_numbers(),
            ase_atoms.get_initial_magnetic_moments()
        )

    return (
        ase_atoms.get_cell(),
        ase_atoms.get_scaled_positions(),
        ase_atoms.get_atomic_numbers()
    )


def get_space_group_symops(atoms: ase.Atoms, use_magnetic_symmetry=False) -> Collection[SymOp]:
    """
    Return the set of symmetry operations that form the space group of a crystal.

    This function uses spglib to extract the space group symmetry operations from an ASE Atoms object,
    and returns them as a collection of `SymOp` objects containing both rotation and translation components.

    Parameters
    ----------
    atoms : ase.Atoms
        An ASE Atoms object representing the crystal structure.

    Returns
    -------
    Collection[SymOp]
        A collection of symmetry operations representing the space group of the structure.
    """
    if use_magnetic_symmetry:
        dataset = spglib.get_magnetic_symmetry(spgcell(atoms, magnetic=True))
    else:
        dataset = spglib.get_symmetry(spgcell(atoms))

    space_group = [SymOp(r, t) for r, t in zip(dataset['rotations'], dataset['translations'])]
    return space_group


def is_group(symops: Collection[SymOp], mod_lattice: bool = False) -> bool:
    """
    Check whether a set of symmetry operations forms a group under composition.

    This function verifies the group axioms:

    * Contains the identity
    * Every element has a (left) inverse in the set
    * Closed under composition

    When `mod_lattice` is True, composition is performed modulo lattice translations.
    This is useful when working with periodic symmetry groups in scaled coordinates.

    Parameters
    ----------
    symops : Collection[SymOp]
        A set or list of symmetry operations to test for group structure.

    mod_lattice : bool, optional
        If True, symmetry operations are composed modulo lattice translations.
        Default is False.

    Returns
    -------
    bool
        True if the set satisfies the group axioms under the specified multiplication rule,
        False otherwise.

    Notes
    -----
    Associativity is assumed due to matrix multiplication being associative.
    """
    if not symops:
        return False

    symop_set = set(symops)

    if len(symop_set) != len(symops):
        return False

    # Identity must be present
    if not any(op.is_identity() for op in symop_set):
        return False

    # Choose the product function based on modulo_lattice
    group_product = symop_product_modulo_lattice if mod_lattice else lambda a, b: a * b

    # Check that every element has an inverse element in the set
    # Note that this has to be done through the products rather than the inverse() member function for SymOps
    # in the case that we calculate modulo lattice translations.
    identity = SymOp.identity()

    for a in symop_set:
        if not any(group_product(a, b).is_identity() for b in symop_set):
            return False

    # Check closure under the group operation
    for a in symop_set:
        for b in symop_set:
            if group_product(a, b) not in symop_set:
                return False

    return True

def is_point_group(symops: Collection['SymOp'], tol: float=magnon.default_numerical_tol) -> bool:
    """
    Check whether a set of symmetry operations forms a point group under composition.

    This function verifies:

    * All symmetry operations are pure rotations (zero translation part).
    * Contains the identity operation.
    * Every element has a (left) inverse in the set.
    * Closed under composition.

    Parameters
    ----------
    symops : Collection[SymOp]
        Set or list of symmetry operations to test.
    tol : float, optional
        Tolerance for checking zero translation. Default is 1e-8.

    Returns
    -------
    bool
        True if the set forms a point group, False otherwise.

    Notes
    -----

    * Point groups must have only rotational parts (no translations).
    * Associativity is assumed due to matrix multiplication being associative.
    """

    if not symops:
        return False

    symop_set = set(symops)

    if len(symop_set) != len(symops):
        return False

    # All operations must have zero translation part
    for op in symop_set:
        if not np.allclose(op.translation, 0, atol=tol):
            return False

    # Identity must be present
    if not any(op.is_identity() for op in symop_set):
        return False

    # Every element has a (left) inverse in the set
    for a in symop_set:
        if not any((a * b).is_identity() for b in symop_set):
            return False

    # Closure under the group operation
    for a in symop_set:
        for b in symop_set:
            if (a * b) not in symop_set:
                return False

    return True


def k_star(k: ArrayLike, symops: Collection[SymOp], tol: float=magnon.default_numerical_tol) -> list[np.ndarray]:
    r"""
    Generate the star of vector k by applying a set of symmetry operations.

    The star of k is defined as:

    .. math::

       \mathrm{star}(\mathbf{k}) = \{ g \cdot \mathbf{k} \mid g \in G \}


    where each symmetry operation g acts on k as a linear transformation, potentially
    followed by a translation.

    Parameters
    ----------
    k : array_like
        A 3-element array representing the input wave vector in scaled coordinates.

    symops : Collection[SymOp]
        A collection of symmetry operations forming a group.

    tol : float, optional
        Tolerance for identifying duplicate vectors in the star. Default is 1e-8.

    Returns
    -------
    list[np.ndarray]
        List of unique wave vectors forming the star of `k`, each as a NumPy array of shape (3,).
    """
    # assert is_point_group(symops)

    star = []
    for s in symops:
        k_new = s @ k
        if not any(np.allclose(k_new, existing, atol=tol) for existing in star):
            star.append(k_new)

    return star


def get_point_group_name(symops: Collection[SymOp]) -> str:
    """
    Return the crystallographic point group name associated with a set of symmetry operations.

    The function extracts the rotational components from the given symmetry operations and uses
    spglib to identify the crystallographic point group to which they belong.

    Parameters
    ----------
    symops : Collection[SymOp]
        A list or set of symmetry operations, typically forming a point group (i.e., all translations removed).

    Returns
    -------
    str
        The international symbol of the crystallographic point group (e.g., "m-3m", "4mm").

    Notes
    -----
    This function assumes the input operations contain only rotation components. If translations are
    present, use `factor_out_translations` or `remove_translations` beforehand.
    """
    rotations = [s.rotation for s in symops]
    return spglib.get_pointgroup(rotations)[0]


def factor_out_translations(symops: Collection[SymOp]) -> SymOp | list[SymOp]:
    """
    Factor out translations from symmetry operations to obtain the associated point group.

    This function takes a list of symmetry operations of the form :math:`(R, \\mathbf{t})`
    and removes the translational components, yielding only the rotational parts :math:`(R, \\mathbf{0})`.
    Duplicate operations (i.e., those with the same rotation matrix) are removed to produce a valid point group.

    Parameters
    ----------
    symops : Collection[SymOp]
        A list or set of symmetry operations (rotation + translation) forming a group.

    Returns
    -------
    list[SymOp]
        The list of unique symmetry operations with zero translations, representing the point group
        associated with the input space group.

    Notes
    -----
    This corresponds to the mathematical quotient of the space group by its translation subgroup:
    :math:`G / T`, where :math:`T` is the group of pure translations.
    """
    point_group = set()
    for s in symops:
        s_rot_only = SymOp(s.rotation, np.zeros(3))
        point_group.add(s_rot_only)

    return list(point_group)

def remove_translations(symops: SymOp | Collection[SymOp]) -> SymOp | list[SymOp]:
    """
    Remove the translational part from one or more symmetry operations.

    This operation replaces each symmetry operation :math:`(R, \\mathbf{t})` with :math:`(R, \\mathbf{0})`,
    effectively projecting the operation into pure rotational space. It is commonly used to convert
    space group operations to their corresponding point group operations.

    Parameters
    ----------
    symops : SymOp or Collection[SymOp]
        A single symmetry operation or a collection of symmetry operations from which to remove translations.

    Returns
    -------
    SymOp | list[SymOp]
        The symmetry operation(s) with translation components set to zero.
        Returns a single SymOp if input is a single object, or a list of SymOps otherwise.

    Notes
    -----
    This function is equivalent to projecting out translations: :math:`(R, \\mathbf{t}) \\rightarrow (R, \\mathbf{0})`.
    If you want to deduplicate the resulting rotations, use `factor_out_translations` instead.
    """
    if isinstance(symops, SymOp):
        return SymOp(symops.rotation, np.zeros(3))

    return [SymOp(s.rotation, np.zeros(3)) for s in symops]


def rotation_axes(symops: Collection[SymOp], tol=magnon.default_numerical_tol) -> list[tuple[np.ndarray, int]]:
    """
    Identify distinct rotation axes from a collection of symmetry operations.

    Each axis is returned as a unit vector along with its associated rotation order. The function excludes
    improper rotations (e.g., mirror planes or inversion) and the identity operation. Axes that are equivalent
    up to sign (i.e., parallel but oppositely directed) are treated as duplicates.

    Parameters
    ----------
    symops : Collection[SymOp]
        A list or set of symmetry operations to analyze.

    tol : float, optional
        Numerical tolerance used to identify equivalent axes (default is 1e-8).

    Returns
    -------
    list[tuple[np.ndarray, int]]
        A list of unique (axis, order) pairs. Each axis is a unit vector (as a NumPy array)
        and the order is the smallest positive integer n such that the rotation raised to the power n equals the identity.
    """
    axes_orders = []
    for s in symops:
        axis, order = s.rotation_axis(), s.rotation_order()
        if axis is None or order == 1:
            continue # skip identity or improper rotations
        if any(np.allclose(axis, a, atol=tol) or np.allclose(axis, -a, atol=tol) for a, o in axes_orders):
            continue
        axes_orders.append((axis, order))

    def sort_key(item):
        axis, order = item
        # Count zeros with tolerance
        num_zeros = np.sum(np.abs(axis) < tol)
        return (-order, -num_zeros)  # Negative for descending sort

    axes_orders_sorted = sorted(axes_orders, key=sort_key)
    return axes_orders_sorted


def site_symmetry_group(space_group: Collection[SymOp], x, cell, delta_equiv_min: float=1e-5) -> list[SymOp]:
    """
    Compute the site symmetry group for a given scaled position `x` in the unit cell.

    Parameters
    ----------
    space_group : Collection[SymOp]
        A list or set of symmetry operations (rotation + translation) forming a space group.

    x : array_like, shape (3,)
        Scaled coordinates of the point for which the site symmetry group is to be determined.

    cell : array_like, shape (3, 3)
        Lattice vectors as rows of a 3×3 matrix defining the unit cell.

    delta_equiv_min : float, optional
        Minimum distance (in Cartesian units) under which two points are considered equivalent.
        Used to identify when a symmetry operation maps the point onto itself modulo lattice translations.
        Default is 1e-5.

    Returns
    -------
    list[SymOp]
        The symmetry operations from the space group that leave `x` fixed (modulo lattice translations).
        These operations may include translations, but they form a point group under multiplication.
    """
    candidate_symops = []
    candidate_squared_distances = []

    x_cartesian = cell @ x
    delta_equiv_min_squared = delta_equiv_min**2

    for si in space_group:
        xs = si @ x
        delta_xs = xs - x
        delta_short_xs = normalize_scaled_coordinate(delta_xs)

        # In Grosse 2002 s_short is written as S_short = S_i X + nearest_int(X - Δ_short_X_S), however, the operation
        # S_i X would produce a vector, so S_short would be a vector not a symmetry operation. I think this is a typo.
        # The text says that S_short should be the symmetry operation that maps X to Δ_short_X_S, so that is what we
        # construct here.
        s_short = si + np.rint(x - delta_short_xs)

        s = s_short
        s_x = s @ x
        s_x_cartesian = cell @ s_x
        x_sx_2 = squared_distance(x_cartesian, s_x_cartesian)

        # find the shortest distance, checking over all adjacent cell translations
        for u in product((-1, 0, 1), repeat=3):
            if u == (0, 0, 0):
                continue

            s_trial = s_short + u
            s_trial_x = s_trial @ x
            s_trial_x_cartesian = cell @ s_trial_x
            delta_trial = squared_distance(x_cartesian, s_trial_x_cartesian)

            if delta_trial < x_sx_2:
                x_sx_2 = delta_trial
                s = s_trial

        # If the initial point and generated equivalent point are close to touching, then add the symmetry
        # operation.
        if x_sx_2 < delta_equiv_min_squared:
            candidate_symops.append(s)
            candidate_squared_distances.append(x_sx_2)

    # sort candidate symops by squared distance
    indices_sorted = sorted(range(len(candidate_squared_distances)), key=lambda i: candidate_squared_distances[i])
    candidate_symops = [candidate_symops[i] for i in indices_sorted]

    site_symmetries = set()
    for s_a in candidate_symops:
        # Ignore duplicate symmetry operations. This is unlikely to ever happen if the input was a proper space group.
        if s_a in site_symmetries:
            continue

        site_symmetries.add(s_a)

        # Check the group is closed without translations modulo lattice
        for s_b in site_symmetries:
            s_ab = s_a * s_b
            if s_ab.is_translation_only(delta_equiv_min):
                site_symmetries.remove(s_a)
                break

    # The site symmetry group with translations factored out must be a point group.
    # assert is_point_group(factor_out_translations(site_symmetries))
    """
    References
    ----------
    R.W. Grosse-Kunstleve et al., "Algorithms for deriving crystallographic space-group information. II. Treatment of
    special positions", Acta Cryst. A **58**, 60–65 (2001). DOI: https://doi.org/10.1107/S0108767301016658

    """

    return list(site_symmetries)


