import magnon
from ase import Atoms
import numpy as np

lattice = np.array([[1, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1]])
positions = np.array([[0, 0, 0], [0.5, 0.5, 0.5]])
magnetic_moments = np.array([[0.0, 0.0, 2.2], [0.0, 0.0, 2.2]])

atoms = Atoms(["Fe", "Fe"], positions=positions, cell=lattice)
atoms.set_initial_magnetic_moments(magnetic_moments)

interactions = magnon.InteractionList([], atoms=atoms)
interactions.insert_interaction(0, 1, np.array([0.5, 0.5, 0.5]), 0.021*np.eye(3))
interactions.insert_interaction(0, 0, np.array([1.0, 0.0, 0.0]), 0.013*np.eye(3))
interactions.insert_interaction(1, 1, np.array([1.0, 0.0, 0.0]), 0.013*np.eye(3))

interactions = interactions.symmetrize(atoms)
atoms, interactions, _ = magnon.build.build_primitive_cell(atoms, interactions)

path = atoms.get_cell().bandpath(path='GHNGPH', npoints=60)

spectrum = magnon.MagnonSpectrum(atoms, interactions)

bstruct = spectrum.get_band_structure(path)
bstruct.plot(emin=0, emax=0.7, filename='BCC_Iron_primitive.png')
