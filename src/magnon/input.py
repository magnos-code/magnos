import re

import ase
import ase.io
import numpy as np

import magnon

def read_site_spin_data(spin_data_filename):
    """
    Read information on the magnitudes and directions of the magnetic moments from a file.

    Parameters
    ----------
        spin_data_filename : str
            The location and name of the file containing the magnitudes and directions of the spin magnetic moments.

        
    Raises
    ----------
        ValueError
            If there are one or more lines which do not match the required format.

        
    Returns
    ----------
        spin_data : numpy.ndarray
            Magnitudes of the spin magnetic moments.
        direction_data : numpy.ndarray
            Directions of the spin magnetic moments.

        
    Notes
    ----------
        The ordering of the lines should correspond to the order of the coordinates in the structure input file.
        Each line consists of the magnitude of the spin magnetic moment followed by the three vector components of its direction.
        These should be separated by spaces.
    """

    spin_data, direction_data = [], []
    with open(spin_data_filename,'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            pattern = r"^([-+]?\d*\.?\d+)\s+([-+]?\d*\.?\d+)\s+([-+]?\d*\.?\d+)\s+([-+]?\d*\.?\d+)"
            match = re.match(pattern, line)

            if match:
                spin_data.append(np.float64(match.group(1)))
                direction_data.append(magnon.linalg.normalised_vector(
                    np.array((match.group(2), match.group(3), match.group(4)), dtype=np.float64)))
                continue

            pattern = r"^([-+]?0*\.?0+)"
            match = re.match(pattern, line)
            if match:
                spin_data.append(np.float64(0.0))
                direction_data.append(np.zeros((3), dtype=np.float64))
                continue

            if not line.isspace():
                raise ValueError(f"Unable to read spin data on line: {line}")
            
    return np.array(spin_data), np.array(direction_data)


def create_interacting_system(cell_file, magmom_file, interaction_file, magmom_scaling=1):
    """
    Read structural, magnetic and exchange coupling data to create an ASE Atoms object and a Magnon.InteractionList object.

    Parameters
    ----------
        cell_file : str
            The location and name of the file containing the cell vectors and coordinates.
        magmom_file : str
            The location and name of the file containing the magnitudes and directions of the magnetic moments.
        interaction_file : str
            The location and name of the file containing the exchange coupling data.
        magmom_scaling : float
            A number by which to multiply all magnetic moment magnitudes.

    Returns
    ----------
        atoms : ase.Atoms object
            The ASE Atoms object describing the structure and magnetic moments.
        interactions : magnon.InteractionList object
            The InteractionList object describing the exchange coupling.

    Notes
    ----------
        The magmom_scaling allows straightforward conversion of the magnitudes if they are provided as e.g. spin quantum numbers
        in the data file (in this case, we would scale by 2 to get the magnetic moment from S).
    """

    atoms = ase.io.read(cell_file)
    magmoms, directions = magnon.input.read_site_spin_data(magmom_file)
    atoms.set_initial_magnetic_moments(magmom_scaling * magmoms[:, np.newaxis] * directions)

    interactions = magnon.interactions.read_interactions(interaction_file)
    interactions = magnon.interactions.InteractionList(interactions, atoms=atoms)

    return atoms, interactions
