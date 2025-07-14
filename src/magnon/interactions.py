import copy
from dataclasses import dataclass

import numpy as np
from ase.utils import reader

from magnon.common import default_distance_tol, find_site_index, lattice_translation_vector

from .symmetry import factor_out_translations, get_space_group_symops, k_star, site_symmetry_group
from .utils import ensure_vector_magnetic_moments


@dataclass
class Interaction:
    """
    A class used to represent a single exchange interaction

    Parameters
    ----------
        i_index : int
            The index of the first site involved in the interaction
        j_index : int
            The index of the second site involved in the interaction
        r_ij_scaled : numpy.ndarray, dtype=float
            The interaction vector from i to j, in units of the lattice vectors
        J_ij : numpy.ndarray, dtype=float
            The exchange coupling tensor
    """

    i_index: int
    j_index: int
    r_ij_scaled: np.ndarray
    J_ij: np.ndarray

    def __iter__(self):
        """
        An iterator which enables unpacking of an Interaction object.
        """
        yield self.i_index
        yield self.j_index
        yield self.r_ij_scaled
        yield self.J_ij


class InteractionList:
    """
    A class used to represent a set of exchange coupling interactions associated with a given structure.
    This class also enables the application of symmetry operations to the exchange coupling dataset.

    Parameters
    ----------
        interaction_data : list
            A list either of Interaction objects or Interaction object initialisation inputs.
        atoms : ase.Atoms object
            An ase.Atoms object describing the structure
        cartesian : bool
            Whether the distance vectors are Cartesian or scaled, defaults to False.

    """

    def __init__(self, interaction_data, atoms, cartesian=False):
        # Ensure that Atoms include magnetic moment vectors (not scalars)
        atoms = ensure_vector_magnetic_moments(atoms)

        self.atoms = atoms.copy()

        self.interaction_list = []
        for i_index, j_index, r_ij, J_ij in interaction_data:
            self.insert_interaction(i_index, j_index, r_ij, J_ij, cartesian)


    def __iter__(self):
        """
        An iterator that yields tuples of interaction data containing i_index, j_index, r_ij, J_ij.
        """
        return iter(self.interaction_list)


    def insert_interaction(self, i_index, j_index, r_ij, J_ij, cartesian=False):
        """
        Appends an Interaction object to the list of interactions.

        Parameters
        ----------
            i_index : int
                The index of the first site involved in the new interaction
            j_index : int
                The index of the second site involved in the new interaction
            r_ij_scaled : numpy.ndarray, dtype=float
                The interaction vector from i to j, in units of the lattice vectors
            J_ij : numpy.ndarray, dtype=float
                The exchange coupling tensor

        """

        if cartesian:
            r_ij_scaled = np.linalg.solve(self.atoms.get_cell().T, r_ij.T).T
        else:
            r_ij_scaled = r_ij

        interaction = Interaction(i_index, j_index, r_ij_scaled, J_ij)
        self.interaction_list.append(interaction)

    def find_interactions(self, i_index, j_index, r_ij_scaled):
        """
        Returns a list of Interaction objects matching the specified indices and vector.

        Parameters
        ----------
            i_index : int
                The index of the first site involved in the interaction(s) to be returned
            j_index : int
                The index of the second site involved in the interaction(s) to be returned
            r_ij_scaled : numpy.ndarray, dtype=float
                The interaction vector from i to j, in units of the lattice vectors, of the interaction(s) to be returned

                
        Returns
        ----------
            matches : list of Interaction objects
                A list of Interaction objects matching the specified indices and vector.

        """

        matches = filter(lambda interaction: interaction.i_index == i_index and interaction.j_index == j_index
                         and np.allclose(interaction.r_ij_scaled, r_ij_scaled), self)
        return list(matches)

    def deduplicate(self):
        """
        Removes any duplicate interactions from the current instance.
        """

        new_interaction_list = []
        coupling_dict = {}
        for interaction in self.interaction_list:
            site_i, site_j, r_ij_scaled, J_ij = interaction

            atom_scaled_positions = self.atoms.get_scaled_positions()

            # Check this is a valid interaction vector
            try:
                u_ij = lattice_translation_vector(r_ij_scaled, atom_scaled_positions[site_i], atom_scaled_positions[site_j])
            except ValueError as e:
                raise RuntimeError(e)

            # The key of site_i, site_j and the lattice translation should be unique
            key = (site_i, site_j, *u_ij)
            if key in coupling_dict:
                if not np.any(np.allclose(coupling_dict[key], J_ij)):
                    raise ValueError(
                        f"Inconsistent value for key {key}: existing value {coupling_dict[key]}, new value {J_ij}")
            else:
                # Only assign if the key does not exist
                coupling_dict[key] = J_ij
                new_interaction_list.append(interaction)

        self.interaction_list = new_interaction_list

    def get_interactions(self, i):
        """
        Returns a list of Interaction objects associated with a specified site index.
        
        Parameters
        ----------
            i : int
                The index of the site to get interactions for.

                
        Returns
        ----------
            A list of Interaction objects in which the first index is i

        """

        return filter(lambda interaction: interaction.i_index == i, self.interaction_list)

    def symmetrize(self,atoms):
        """
        Apply the symmetry of the cell to the couplings

        
        Parameters
        ----------
            atoms : ASE.Atoms object
                The ASE Atoms object of the cell whose symmetry to use.

                
        Returns
        ----------
            InteractionList object containing the symmetrized exchange coupling

                
        See Also
        ----------
            magnon.interactions.symmetrize_couplings

        """
        interaction_data = symmetrize_couplings(self, atoms)

        return InteractionList(interaction_data, atoms, cartesian=False)

    def update_indices(self,index_map):
        """
        Update the site indices based on a provided remapping
                
        Parameters
        ----------
            index_map : list
                A list in which the entry at each position is the new site index.

        Raises
        ----------
            ValueError
                If one or more of the updated indices is invalid.

        Returns
        ----------
            InteractionList object containing the updated site indices.
   
        Notes
        ----------
            This can be used to update the site indices when the nonmagnetic sites are removed - see `magnon.MagnonSpectrum`.
            The mapping in this case contains -1 for the nonmagnetic sites. If -1 is detected in the new couplings, a ValueError is raised.

        """
        new_interaction_list = []
        for interaction in self.interaction_list:
            site_i, site_j, r_ij_scaled, J_ij = interaction
            new_i_index = index_map[site_i]
            new_j_index = index_map[site_j]
            if new_i_index < 0:
                raise ValueError(f"The new index for {site_i} is invalid")
            if new_j_index < 0:
                raise ValueError(f"The new index for {site_j} is invalid")
            new_interaction_list.append(Interaction(new_i_index, new_j_index, r_ij_scaled, J_ij))

        return InteractionList(new_interaction_list, self.atoms)

def apply_bond_reversal_symmetry(interaction_list: InteractionList):
    """
    Applies the bond reversal symmetry due to bidirectionality of the exchange coupling.
    Returns the new couplings generated under this symmetry which were not included in the original interaction data.

    Parameters
    ----------
        interaction_list : InteractionList object
            The exchange coupling to which the bond reversal symmetry will be applied.

    Raises
    ----------
        RuntimeError
            If the symmetry results in inconsistent interactions.
        
    Returns
    ----------
        new_interaction_list : InteractionList object
            The new interactions generated under this symmetry which were not included in the original interaction data.
        
    Notes
    ----------
        Applies the intrinsic symmetry

        .. math:: J(r_ij) = (J(r_ji))^T

        Any symmetric couplings which are absent in the input are returned. If any existing couplings violate the intrinsic symmetry then a RuntimeError is raised.
        This will happen in the following cases:
        - More than 1 interaction is found that has the same :math:`i,j` -> :math:`j,i` and :math:`r_ij` -> :math:`r_ji` symmetry. This means
        the input data is incorrect or incompatible with our code (for example assuming that identical interactions
        will be summed).
        - An interaction is found with the :math:`i,j` -> :math:`j,i` and :math:`r_ij` -> :math:`r_ji` symmetry but the coupling matrix does not have
        the correct :math:`J(r_ij) = (J(r_ji))^T` symmetry.

    """

    new_interaction_list = copy.deepcopy(interaction_list)

    for interaction in interaction_list:

        # This is what the intrinsic symmetric interaction should look like
        sym_interaction = Interaction(
            interaction.j_index, interaction.i_index, -interaction.r_ij_scaled, np.transpose(interaction.J_ij))

        matches = new_interaction_list.find_interactions(sym_interaction.i_index, sym_interaction.j_index, sym_interaction.r_ij_scaled)

        # if no symmetric interaction was found for the site indices and the interaction vector then insert it as a new
        # interaction
        if len(matches) == 0:
            new_interaction_list.insert_interaction(sym_interaction.i_index, sym_interaction.j_index, sym_interaction.r_ij_scaled, sym_interaction.J_ij)
            continue
        elif len(matches) == 1:
            # if a symmetric interaction was found, check the exchange matrix obeys the required symmetry
            if not np.allclose(matches[0].J_ij, sym_interaction.J_ij):
                raise RuntimeError(
                    f"Intrinsic symmetry J(r_ij) = (J(r_ji))^T is broken because the symmetric interaction did not "
                    f"obey the transpose symmetry.\n"
                    f"Initial interaction is   {interaction.i_index} {interaction.j_index} {interaction.r_ij_scaled}, {interaction.J_ij.flatten()}\n"
                    f"Symmetric interaction is {sym_interaction.i_index} {sym_interaction.j_index} {sym_interaction.r_ij_scaled}, {sym_interaction.J_ij.flatten()}")
        else:
            # more than one symmetric interaction was found, this should not be possible,
            # there is probably something wrong with the input data
            raise RuntimeError(
                "Intrinsic symmetry check J(r_ij) = (J(r_ji))^T found more than one symmetric interaction")

    return new_interaction_list


def symmetrize_couplings(interaction_list, atoms, distance_tol: float=default_distance_tol):
    """
    Generates the full exchange interaction data under the symmetry of a unit cell

    Parameters
    ----------
        interaction_list : magnon.InteractionList object
            The InteractionList object containing the exchange coupling to be symmetrized
        atoms : ASE.Atoms object
            The ASE Atoms object of the cell associated with the interactions
        distance_tol : float, optional
            The precision used in determining equivalent coordinates, defaults to 1e-6

        
    Raises
    ----------
        RuntimeError
            If the image of an atom under a symmetry operation does not correspond to another atom's location

        
    Returns
    ----------
        new_interaction_list : magnon.InteractionList object
            The fully-symmetrized exchange coupling interactions
    """


    #Get symmetries of input cell
    space_group_symops = get_space_group_symops(atoms)

    site_symmetries = []
    for site in atoms.get_scaled_positions():
        site_symmetries.append(factor_out_translations(site_symmetry_group(space_group_symops, site, atoms.get_cell())))

    interaction_list = apply_bond_reversal_symmetry(interaction_list)

    new_interaction_list = copy.deepcopy(interaction_list)
    scaled_positions = atoms.get_scaled_positions()
    for (site_i, site_j, r_ij, J_ij) in interaction_list:
        r_i, r_j = scaled_positions[site_i], scaled_positions[site_j]

        for new_r_ij in k_star(r_ij, site_symmetries[site_i]):
            new_r_j = new_r_ij + r_i
            new_site_j = find_site_index(new_r_j, scaled_positions, distance_tol=distance_tol)

            if new_site_j is None:
                raise RuntimeError("failed to find symmetric site_i")

            new_interaction_list.insert_interaction(site_i, new_site_j, new_r_ij, J_ij)


    new_interaction_list.deduplicate()

    return new_interaction_list


@reader
def read_interactions(filename, atom_base_index=1):
    """
    Reads in interaction data from a file.

    Parameters
    ----------
        filename : str
            The location of the file containing the interaction data.
        atom_base_index : int, optional
            The convention for the index of the first atom, defaults to 1
        
    Raises
    ----------
        ValueError
            If the number of entries does not enable the construction of a scalar or 3x3 exchange coupling
        
    Returns
    ----------
        interaction_data : list
            A list of data needed to create an Interaction object for each exchange coupling interaction in the file.
        
    Notes
    ----------
        Each line describes an interaction. The first two values are the atom indices
        in the unit cell, the next three give the vector between the atoms and the remainder give the coupling strength, either
        as a tensor or a scalar.
    """

    fd = filename

    interaction_data = []

    line = fd.readline()
    while line:
        if line.startswith('#') or line.startswith('//') or line.startswith('%') or line.isspace():
            line = fd.readline()
            continue

        line_parts = line.split()

        i_index = int(line_parts[0]) - atom_base_index
        j_index = int(line_parts[1]) - atom_base_index

        r_ij = np.array([float(line_parts[2]), float(line_parts[3]), float(line_parts[4])])

        if len(line_parts[5:]) == 1:
            # scalar exchange
            J_ij = float(line_parts[5]) * np.identity(3)
        elif len(line_parts[5:]) == 9:
            # tensor exchange
            tensor = np.array([float(x) for x in line_parts[5:]])
            J_ij = np.reshape(tensor, (3, 3))
        else:
            raise ValueError(f'Invalid exchange line: {line}')

        interaction_data.append((i_index, j_index, r_ij, J_ij))

        line = fd.readline()

    return interaction_data
