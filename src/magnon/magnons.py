import multiprocessing

import ase.dft
import ase.spectrum.band_structure
import ase.spectrum.dosdata
import numpy as np
from numpy.typing import ArrayLike

import magnon


class MagnonSpectrum:
    """
    A class to build and diagonalize the Linear Spin Wave Hamiltonian from exchange interactions.

    Parameters
    ----------
        atoms : ase.Atoms object
            The Atoms object of the system to which the magnon spectrum pertains.
        interactions : magnon.interactions.InteractionList object
            The InteractionList object of the exchange interactions to use for the spectrum calculation.
        num_threads : int
            The number of threads to use for parallel calculations.

    Other Parameters
    ----------------
        mag_order_vect : numpy.ndarray, dtype: float
            Wavevector for Q-commensurate spin vector ordering.

        mag_order_axis : numpy.ndarray, dtype: float
            The rotation axis around which the Q-commensurate spin vector ordering rotations are to occur

        ham_prefactor : float
           Choice of Hamiltonian prefactor convention.

        spins_are_unit : bool
            Whether to use the unit spin model, or dimensionful spins.

        supercell : array-like
            Size of the calculation cell with respect to the primitive cell

        name : str
            The name of the calculation

    Attributes
    ----------
        u : numpy.ndarray, dtype: float
            The u vectors for each site. See :func:`magnon.magnons.MagnonSpectrum._generate_uv_vectors`.

        v : numpy.ndarray, dtype: float
            The v vectors for each site. See :func:`magnon.magnons.MagnonSpectrum._generate_uv_vectors`.

        S : numpy.ndarray, dtype: float
            Spin quantum numbers.

        usm_factors : numpy.ndarray, dtype: float
           Unit spin model conversion factors.

    See Also
    ----------

    This class mirrors the structure of ASE.phonons.Phonons.

    """


    def __init__(self, atoms: ase.Atoms, interactions: magnon.interactions.InteractionList, calc=None,
                 supercell: ArrayLike = (1, 1, 1), name: str | None = None, num_threads: int = 1, **kwargs):
        self._num_threads = num_threads

        if 'name' not in kwargs:
            kwargs['name'] = "magnon"

        # Check atoms has magmoms and that it's vectorial.
        atoms = magnon.utils.ensure_vector_magnetic_moments(atoms)

        # indices of atoms which have non-zero magnetic moments
        if atoms.get_initial_magnetic_moments().ndim == 2:
            mag_indices = np.where(atoms.get_initial_magnetic_moments() != np.zeros(3))[0]
        elif atoms.get_initial_magnetic_moments().ndim == 1:
            mag_indices = np.where(atoms.get_initial_magnetic_moments() != 0)[0]
        else:
            raise ValueError("Magnetic moments must be scalar or vector quantities.")

        if len(mag_indices) == 0:
            raise ValueError("No magnetic moments in the 'atoms' input")

        # Map indices from all sites to magnetic sites only
        map_to_mag_indices = -1*np.ones((len(atoms)))
        for mag,nonmag in enumerate(mag_indices):
            map_to_mag_indices[nonmag] = mag

        interactions = interactions.update_indices(mag_indices)
        mag_atoms = atoms[mag_indices]

        self.interactions = interactions
        self.atoms = mag_atoms

        self.cell_vectors = None
        self.interaction_matrix = None
        self._build_interaction_matrix()

        self.calc = calc


        self.name = name
        self.supercell = supercell

        self.u, self.v = MagnonSpectrum._generate_uv_vectors(mag_atoms.get_initial_magnetic_moments())

        if 'mag_order_vect' in kwargs:
            self.mag_order_vect = kwargs['mag_order_vect']
        else:
            self.mag_order_vect = np.zeros(3)

        if 'mag_order_axis' in kwargs:
            self.mag_order_axis = kwargs['mag_order_axis']
        else:
            self.mag_order_axis = np.array([0, 1, 0])

        if 'ham_prefactor' in kwargs:
            self.ham_prefactor = kwargs['ham_prefactor']
        else:
            self.ham_prefactor = 2.0

        if 'spins_are_unit' in kwargs:
            self.spins_are_unit = kwargs['spins_are_unit']
        else:
            self.spins_are_unit = True

        self.spin = 0.5*np.linalg.norm(mag_atoms.get_initial_magnetic_moments(), axis=1)
        self.S, self.usm_factors = magnon.lattice.unit_spin_model_factors(self.spin, self.spins_are_unit)


    def _build_interaction_matrix(self):
        mag_atoms = self.atoms
        nmag_atoms = len(mag_atoms)

        cartesian_shift_vectors = []
        for _, _, r_ij_scaled, _ in self.interactions:
            cartesian_shift_vectors.append(self.atoms.get_cell().T @ r_ij_scaled)
        cartesian_shift_vectors = np.array(cartesian_shift_vectors)

        cell_vectors, cell_vector_lookup = np.unique(cartesian_shift_vectors, axis=0, return_inverse=True)

        interaction_matrix = np.zeros((cell_vectors.shape[0], nmag_atoms, nmag_atoms, 3, 3), dtype=np.float64)
        for n, (i, j, _, J_ij) in enumerate(self.interactions):
            shift_index = cell_vector_lookup[n]
            interaction_matrix[shift_index, i, j, :, :] = J_ij

        self.cell_vectors = np.array(cell_vectors)
        self.interaction_matrix = np.array(interaction_matrix)

    @staticmethod
    def _generate_uv_vectors(spin_vectors: np.ndarray) -> (np.ndarray, np.ndarray):
        """
        Compute u and v vectors from spin orientation vectors for LSWT calculation
        """

        spin_unit_vectors = np.einsum("ij,i->ij", spin_vectors,  1/np.linalg.norm(np.atleast_2d(spin_vectors), axis=1))
        u_vectors, v_vectors = [], []
        z = np.array([0, 0, 1.0])
        u_z = np.array([1.0, 1.0j, 0])
        u_minus_z = np.array([1, -1.0j, 0])
        for i in range(spin_unit_vectors.shape[0]):
            dot_z = np.dot(z, spin_unit_vectors[i].astype(np.float64))
            if np.isclose(dot_z, 1):
                u_vectors.append(u_z)
            elif np.isclose(dot_z, -1):
                u_vectors.append(u_minus_z)
            else:
                # Get rotation matrix
                axis = np.cross(z, spin_unit_vectors[i])
                axis /= np.dot(axis, axis)
                angle = np.arccos(dot_z)
                rot, _ = magnon.linalg.rotation_matrix_pair(axis, angle)
                u_vectors.append(np.einsum("ij,j->i", rot, u_z))
            v_vectors.append(spin_unit_vectors[i])

        return np.array(u_vectors, dtype=np.complex128), np.array(v_vectors, dtype=np.complex128)

    def calculate_eigen_energy(self, k: ArrayLike) -> np.ndarray:
        """
        Compute the eigenvalues of the Hamiltonian at wavevector k. The Hamiltonian in the rotated frame approach [1]_ is built and
        diagonalized using the treatment for generalized quadratic bosonic Hamiltonians [2]_.
        
        Parameters
        ----------
            k : ndarray, shape (3,)
                Wavevector in reciprocal space.
                
        Returns
        ----------
            hamiltonian_eigs : ndarray
                Eigenvalues of the Bogoliubov-transformed Hamiltonian.

        References
        ----------

        .. [1] S. Toth and B. Lake, J. Phys. Condens. Matter, 27, 166002 (2015)
        .. [2] J.H.P Colpa, Physica, 93A, 327-353 (1978)

        """
        mag_atoms = self.atoms
        nmag_atoms = len(mag_atoms)

        # Prepare variables outside the X loop
        u = self.u
        ubar = np.conjugate(u)
        sqrt_s = np.sqrt(self.S)
        sqrt_usm = np.sqrt(self.usm_factors)
        exp_k_r = np.exp(1j * np.dot(mag_atoms.get_positions(), k))

        # Allocate term arrays
        term_1 = np.zeros((nmag_atoms, nmag_atoms), dtype=np.complex128)
        term_2 = np.zeros((nmag_atoms, nmag_atoms), dtype=np.complex128)
        term_3 = np.zeros((nmag_atoms, nmag_atoms), dtype=np.complex128)
        term_4 = np.zeros((nmag_atoms, nmag_atoms), dtype=np.complex128)
        term_z = np.zeros((nmag_atoms, nmag_atoms), dtype=np.complex128)

        inv_cell = np.linalg.inv(self.atoms.get_cell())

        for x_index, x in enumerate(self.cell_vectors):
            # Precompute rotation matrix and exponentials for this x
            x_scaled = inv_cell @ x
            exp_k_x = np.exp(1j * np.dot(k, x))
            Q_cell_ang = 2 * np.pi * np.dot(self.mag_order_vect, x_scaled)
            Q_cell = magnon.linalg.rotation_matrix_pair(self.mag_order_axis, Q_cell_ang)[0]

            for b in range(nmag_atoms):
                for b_prime in range(nmag_atoms):
                    J_x_bb = self.interaction_matrix[x_index, b, b_prime] * sqrt_usm[b] * sqrt_usm[b_prime]
                    J_x_bb_rot = J_x_bb @ Q_cell

                    prefactor = -0.25 * sqrt_s[b] * sqrt_s[b_prime] * exp_k_x * exp_k_r[b_prime] * np.conjugate(
                        exp_k_r[b])

                    term_1[b, b_prime] += prefactor * np.einsum('i,ij,j', u[b], J_x_bb_rot, ubar[b_prime])
                    term_2[b, b_prime] += prefactor * np.einsum('i,ij,j', u[b], J_x_bb_rot, u[b_prime])
                    term_3[b, b_prime] += prefactor * np.einsum('i,ij,j', ubar[b], J_x_bb_rot, ubar[b_prime])
                    term_4[b, b_prime] += prefactor * np.einsum('i,ij,j', ubar[b], J_x_bb_rot, u[b_prime])

                    if b == b_prime:
                        # Only needed on diagonal
                        z_sum = 0
                        for b_dblprime in range(nmag_atoms):
                            J_x_bb = self.interaction_matrix[x_index, b_dblprime, b_prime] * sqrt_usm[b_prime] * sqrt_usm[b]
                            J_x_bb_rot = J_x_bb @ Q_cell
                            z_sum += 0.5 * self.S[b_dblprime] * np.einsum(
                                'i,ij,j', self.v[b_dblprime], J_x_bb_rot, self.v[b_prime]
                            )
                        term_z[b, b_prime] += z_sum

        # Add z term to correct blocks
        term_1 += term_z
        term_4 += term_z

        # Hamiltonian assembly
        eps = 1e-4
        eps_block = eps * np.identity(nmag_atoms, dtype=np.complex128)
        small_ham = self.ham_prefactor * np.block([
            [term_1 + eps_block, term_2],
            [term_3, term_4 + eps_block]
        ])

        paraidentity = magnon.linalg.paraidentity(nmag_atoms)

        # Bogoliubov transform
        try:
            ham_chol_lower = np.linalg.cholesky(small_ham)
            ham_chol_upper = np.conjugate(ham_chol_lower.T)
        except np.linalg.LinAlgError:
            hamiltonian_eigs, _ = np.linalg.eig(small_ham)
            print("Eigenvalues of small ham are:")
            print(hamiltonian_eigs)
            raise Exception("[ERROR] Cholesky decomposition failed - check that the Hamiltonian is positive definite.")

        upper_para_lower = ham_chol_upper @ paraidentity @ ham_chol_lower
        upl_eigs, upl_eigvects = np.linalg.eig(upper_para_lower)
        upl_argsort = np.flip(np.argsort(upl_eigs))

        U = upl_eigvects[:, upl_argsort]
        L = np.diag(upl_eigs[upl_argsort])

        E = paraidentity @ L
        hamiltonian_eigs = np.diagonal(E) - eps

        return hamiltonian_eigs

    def get_band_structure(self, path: ase.dft.kpoints.BandPath) -> ase.spectrum.band_structure.BandStructure:
        """
        Gets the magnon bandstructure along a path given by an ASE BandPath object. This is a wrapper for `band_structure()`
        which enables handling of k-space paths and bandstructures using ASE objects.

        Parameters
        ----------
            path : ASE.dft.kpoints.BandPath object
                The path along which to compute the bandstructure.

        Returns
        -------
            ase.spectrum.band_structure.BandStructure
                An ASE BandStructure object containing the magnon bandstructure.

        """
        energy_kl = self.band_structure(path.cartesian_kpts())
        bs = ase.spectrum.band_structure.BandStructure(path, energies=energy_kl)
        return bs

    def band_structure(self, path_kc: np.ndarray) -> np.ndarray:
        """
        Computes the bands energies for a set of reciprocal space points.

        Parameters
        ----------
            path_kc : array-like
                The reciprocal space path along which to compute the bandstructure, given in Cartesian coordinates.

        Returns
        -------
            numpy.ndarray
                The band energies along the specified path.

        """
        assert self.interactions is not None

        energy_kl = []
        if self._num_threads == 1:
            # Doing pool processing for 1 thread adds too much overhead, so call the function directly.
            energy_kl = [_single_q_energy_worker((self, q_c)) for q_c in path_kc]
        else:
            with multiprocessing.Pool(processes=self._num_threads) as pool:
                energy_kl = list(pool.imap(_single_q_energy_worker, [(self, q_c) for q_c in path_kc]))

        return np.array([energy_kl])

    def get_dos(self, kpts: ArrayLike=(10, 10, 10)) -> ase.spectrum.dosdata.RawDOSData:
        """
        Computes the magnon density of states.

        Parameters
        ----------
            kpts : array-like, dtype=int
                The dimensions of the Monkhorst-Pack grid to sample.

        Returns
        -------
            ase.spectrum.dosdata.RawDOSData object

        """
        kpts_kc = ase.dft.kpoints.monkhorst_pack(kpts)
        w_w = self.band_structure(kpts_kc).ravel()
        dos = ase.spectrum.dosdata.RawDOSData(w_w, np.ones_like(w_w))
        return dos

def _single_q_energy_worker(args: (MagnonSpectrum, np.ndarray)) -> np.ndarray:
        """
        Compute and return sorted real magnon energies at a k-point.

        This is a multiprocessing-friendly worker function for use with MagnonSpectrum.
        It must be defined at module level to be pickleable.

        
        Parameters
        ----------
            args : tuple
                Tuple of (instance, q_c) where:
                    instance : MagnonSpectrum
                        Instance of the MagnonSpectrum class.
                    q_c : array_like
                        k-point in scaled Cartesian coordinates.

                
        Returns
        ----------
            energies : ndarray
                Sorted real eigenvalues of the magnon Hamiltonian at the given k-point.

                
        Notes
        ----------
            The k-point is converted to cartesian space (multiplied by 2*pi).
            Only the real parts of the eigenvalues are returned.
        """
        instance, q_c = args
        k_cartesian = 2 * np.pi * np.array(q_c)
        energy = instance.calculate_eigen_energy(k_cartesian)
        energy.sort()
        return energy.real


