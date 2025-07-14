import magnon
from ase import Atoms
import numpy as np

lattice = np.array([[1,0,0],
                    [0,1,0],
                    [0,0,1]])
positions = np.array([[0,0,0],[0.5,0.5,0.5]])
magnetic_moments = np.array([[0,0,1], [0,0,1]])

atoms = Atoms(["Fe", "Fe"], positions=positions, cell=lattice)
atoms.set_initial_magnetic_moments(magnetic_moments)

interactions = magnon.InteractionList([], atoms=atoms)
interactions.insert_interaction(0, 1, np.array([0.5, 0.5, 0.5]), np.eye(3))

interactions = magnon.interactions.apply_bond_reversal_symmetry(interactions)

print("Bond reversal symmetry")
for i,j,r_ij,J_ij in interactions:
    print(f"i ={i:2d}, j ={j:2d}, r_ij = ({r_ij[0]:4.1f}, {r_ij[1]:4.1f}, {r_ij[2]:4.1f}), J_ij ={J_ij[0,0]:5.2f} eV")
print()

interactions = interactions.symmetrize(atoms)

print("Symmetry operations")
for i,j,r_ij,J_ij in interactions:
    print(f"i ={i:2d}, j ={j:2d}, r_ij = ({r_ij[0]:4.1f}, {r_ij[1]:4.1f}, {r_ij[2]:4.1f}), J_ij ={J_ij[0,0]:5.2f} eV")