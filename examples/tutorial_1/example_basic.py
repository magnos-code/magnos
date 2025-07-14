import magnon
from ase import Atoms
import numpy as np

lattice = np.array([[1,0,0],
                    [0,1,0],
                    [0,0,1]])
positions = np.array([[0.5,0.5,0.5]])
magnetic_moments = np.array([[0,0,1]])

atoms = Atoms("Fe", positions=positions, cell=lattice)
atoms.set_initial_magnetic_moments(magnetic_moments)

interactions = magnon.InteractionList([], atoms=atoms)
interactions.insert_interaction(0, 0, np.array([-1, 0 ,0]), np.eye(3))
interactions.insert_interaction(0, 0, np.array([ 1, 0, 0]), np.eye(3))

special_points = {
    'G' : [0, 0, 0],
    'X' : [0.5, 0, 0],
}

path = atoms.get_cell().bandpath(path='XGX', npoints=100, special_points=special_points)

spectrum = magnon.MagnonSpectrum(atoms, interactions)

bstruct = spectrum.get_band_structure(path)
bstruct.plot(emin=0, emax=20, filename='0basic_bands.png')
